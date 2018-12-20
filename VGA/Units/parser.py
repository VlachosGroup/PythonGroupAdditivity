import re
from ..Error import UnitsParseError
from .db import units_db


__all__ = ['eval_expr']


class UnitsParser(object):
    # Regular expression for tokenizing units declarations.
    # TODO: Support evaluation of numbers with e/E exponents.
    tokenize_re = re.compile(r'-?[.\d]+|[a-zA-Z]+|.')

    def __init__(self, expr, debug=False):
        self.tokens = [tok for tok in re.findall(self.tokenize_re, expr)
                       if not tok.isspace()]
        # print 'tokens:', self.tokens
        self.idx = 0
        self.depth = 0
        self.debug = debug

    def isnumber(self, what):
        try:
            float(what)
            return True
        except ValueError:
            return False

    # These two are used to output debugging info.
    def enter(self, what):
        if self.debug:
            print('%s%r {' % (' '*(self.depth*4), what))
            self.depth += 1

    def leave(self, what):
        if self.debug:
            self.depth -= 1
            print('%s} %r' % (' '*(self.depth*4), what))

    def peek(self):
        if self.idx == len(self.tokens):
            return None
        return self.tokens[self.idx]

    def take(self):
        val = self.peek()
        if val is not None:
            self.idx += 1
        return val

    def fork(self):
        return self.idx

    def join(self, old_idx):
        self.idx = old_idx

    def parse(self):
        self.enter('[root]')
        result = None
        try:
            result = self.parse_expr()
            next = self.peek()
            if next is not None:
                result = None
                raise UnitsParseError("Unexpected token '%s'\n" % next)
        finally:
            self.leave(result)
        return result

    def parse_expr(self):
        self.enter('expr')
        result = None
        try:
            result = ('expr', self.parse_factor())
            next = self.peek()
            while True:
                if next is None:
                    # result = None
                    break
                elif next in '*/':
                    self.take()
                    result = ('expr', result, next, self.parse_factor())
                    next = self.peek()
                else:
                    old_state = self.fork()
                    try:
                        result = ('expr', result, '*', self.parse_factor())
                        next = self.peek()
                    except UnitsParseError:
                        self.join(old_state)
                        break
        finally:
            self.leave(result)
        return result

    def parse_factor(self):
        self.enter('factor')
        result = None
        try:
            left = self.parse_base()
            next = self.peek()
            if next == '^':
                self.take()
                result = ('factor', left, '^', self.parse_number())
            else:
                result = ('factor', left)

        finally:
            self.leave(result)
        return result

    def parse_base(self):
        self.enter('base')
        result = None
        try:
            next = self.peek()
            if next is None:
                result = None
                raise UnitsParseError(
                    "Expected name or number, got end of input instead")
            elif next == '(':
                self.take()
                expr = self.parse_expr()
                if self.take() != ')':
                    raise UnitsParseError(
                        "Expected closing parenthesis, got '%s' instead"
                        % next)
                result = ('base', expr)
            elif self.isnumber(next):
                result = ('base', self.parse_number())
            else:
                result = ('base', self.parse_name())

        finally:
            self.leave(result)
        return result

    def parse_name(self):
        self.enter('name')
        result = None
        try:
            next = self.take()
            if next is None:
                raise UnitsParseError(
                    "Expected name instead of end of input")
            if not next.isalpha():
                raise UnitsParseError(
                    "Expected name, got '%s' instead" % next)
            result = ('name', next)

        finally:
            self.leave(result)
        return result

    def parse_number(self):
        self.enter('number')
        result = None
        try:
            next = self.take()
            if next == '(':
                next = self.take()
                closing = self.take()
                if closing != ')':
                    raise UnitsParseError("Expected closing parenthesis, "
                                          "got '%s' instead" % closing)
            if next is None:
                raise UnitsParseError(
                    "Expected number instead of end of input")
            if self.isnumber(next):
                if '.' in next:
                    result = ('number', float(next))
                else:
                    result = ('number', int(next))
            else:
                raise UnitsParseError(
                    "Expected number, got '%s' instead" % next)

        finally:
            self.leave(result)
        return result


def parse(expr, *args, **kwargs):
    """Parse units expression (expr) and return abstract syntax tree."""
    return UnitsParser(expr, *args, **kwargs).parse()


def eval_subtree(tree):
    """
    Evaluate subtree of abstract syntax tree (obtained from
    :func:`parser.parse`).
    """
    if tree[0] == 'expr':
        # Subtree represents an *expression*.
        if len(tree) > 2:
            if tree[2] == '*':
                # Multiply two *factors*.
                return eval_subtree(tree[1])*eval_subtree(tree[3])
            elif tree[2] == '/':
                # Divide two *factors*.
                return eval_subtree(tree[1])/eval_subtree(tree[3])
        else:
            return eval_subtree(tree[1])

    elif tree[0] == 'factor':
        # Subtree represents a *factor* in an expression.
        if len(tree) > 2:
            if tree[2] == '^':
                return eval_subtree(tree[1])**eval_subtree(tree[3])
        else:
            return eval_subtree(tree[1])

    elif tree[0] == 'base':
        # Subtree represents the *base* of an exponentiation.
        return eval_subtree(tree[1])

    elif tree[0] == 'name':
        # Subtree represents *name* of a unit; lookup unit and return.
        return units_db.lookup(tree[1])

    elif tree[0] == 'number':
        # Subtree represents a *number*.
        return tree[1]

    # This shouldn't happen...
    assert False


def eval_expr(expr):
    """
    Evaluate string expression of physical quantity or units.
    """
    ast = parse(expr)
    return eval_subtree(ast)
