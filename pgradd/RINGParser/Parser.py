from ..Error import RINGSyntaxError
from operator import eq

"""
Parses the RING input string into abstract syntax tree as defined by
Grammar.py.

Example:
>>>from pgradd.RINGParser.parser import parse
>>>parse('\
>>>fragment C3Chain{\
>>>   C labeled c1\
>>>   C labeled c2 single bond to c1\
>>>   C labeled c3 single bond to c2\
>>>   }')
Out:
[RINGToken('RINGInput'),
 [RINGToken('Fragment'),
  [RINGToken('FragmentName'), 'C3Chain'],
  [RINGToken('MolQuery'),
   [RINGToken('Atom'),
    [RINGToken('AtomType'), [RINGToken('Symbols'), 'C']],
    [RINGToken('AtomLabel'), 'c1']],
   [RINGToken('AtomChain'),
    [RINGToken('BondedAtom'),
     [RINGToken('AtomType'), [RINGToken('Symbols'), 'C']],
     [RINGToken('AtomLabel'), 'c2'],
     [RINGToken('BondType'), 'single'],
     [RINGToken('AtomLabel'), 'c1']],
    [RINGToken('AtomChain'),
     [RINGToken('BondedAtom'),
      [RINGToken('AtomType'), [RINGToken('Symbols'), 'C']],
      [RINGToken('AtomLabel'), 'c3'],
      [RINGToken('BondType'), 'single'],
      [RINGToken('AtomLabel'), 'c2']]]]]]]
"""


class RINGToken(object):
    # Each node is converted to RINGToken
    def __init__(self, name):
        self.name = name

# TODO: Resolve unknown name 'cmp'

    def __cmp__(self, other):
        if isinstance(other, RINGToken):
            return eq(self.name, other.name)
        else:
            return eq(self.name, other)

    def __str__(self):
        return self.name

    def __repr__(self):
        return "RINGToken('%s')" % self.name


class Parser(object):
    # Parser is a parent class for various type of syntax
    # syntax classes:
    # EOS               : end of the line identifier
    # Digit             : 1 digit class
    # Number            : mutiple digit class
    # String            : string class. (when to stop interpret is defined
    #                     by string_stop
    # Literal           : Has to match the specified symbol to be identified
    #                     as a class
    # Filler            : Looks for the match, but the information is
    #                     not stored.
    # DeprecatedLiteral : similar to Literal, but warning is raised.
    # Optional          : Indiciates that child class of this class is
    #                     not necessary
    # All               : Interpreter goes through all the child classes
    #                     of this class
    # Either            : At least one of the child class has to be in
    #                     the stream
    # ZeroOrMore        : Similar to Optional but can take multiple input
    # Literals          : Has to match one of the strings in this class.
    #                     Equal to Either(Literal('something'),
    #                     Literal('Something'),...)
    #
    # Notes:
    # The code compares the stream with given syntax tree in the grammar. so
    # ParseState.parse() stacks with the number of child class. When the errors
    # occur, the parser undo each parse loop by using 'try'. Eventually,
    # The error is given at the user's calling function.
    # For classes like Either and Optional, the algorithm still uses errors
    # when there is no match, but the error is suppressed using 'with'
    # function, A lot of valuable coding is in this work.
    def set_name(self, name):
        self.name = name


class EOS(Parser):
    def __call__(self, stream, output):
        if stream.peek() != '':
            stream.error('<end of string>')

    def __str__(self):
        return '<end of string>'

    def __repr__(self):
        return 'EOS()'


class Digit(Parser):
    def __init__(self, n=1):
        self.n = n

    def __call__(self, stream, output):
        out = ''
        for i in range(self.n):
            c = stream.peek()
            if not c.isdigit():
                stream.error('<digit>')
            out += stream.take()
        output.append(int(out))

    def __str__(self):
        if self.n == 1:
            return '<digit>'
        else:
            return '<%d digits>' % self.n

    def __repr__(self):
        if self.n == 1:
            return 'Digit()'
        else:
            return 'Digit(%r)' % self.n


class Number(Parser):
    def __init__(self):
        pass

    def __call__(self, stream, output):
        out = stream.peek()
        if not out.isdigit():
            stream.error('<number>')
        stream.take()
        while stream.peek().isdigit():
            out += stream.take()
        output.append(int(out))

    def __str__(self):
        return '<number>'

    def __repr__(self):
        return 'Number()'


string_okay = ['_']


class String(Parser):
    # Read continuous input of alphabet, digit, and any other character in
    # string_okay
    def __init__(self):
        pass

    def __call__(self, stream, output):
        if not (stream.peek().isalpha() or
                stream.peek().isdigit() or
                stream.peek() in string_okay):
                stream.error('<string>')
        nn = 2

        while (stream.peek(nn)[-1].isalpha() or
               stream.peek(nn)[-1].isdigit() or
               stream.peek(nn)[-1] in string_okay):
            nn += 1
        out = stream.take(n=nn-1)
        output.append(out)

    def __str__(self):
        return '<string>'

    def __repr__(self):
        return 'String()'


class Literal(Parser):
    def __init__(self, tok, no_error=False):
        self.tok = tok
        self.no_error = no_error

    def __call__(self, stream, output):
        if stream.peek(len(self.tok)) != self.tok:
            if self.no_error:
                # Suppress error message.  (It's handled by something further
                # up the stack).
                stream.error(None)
            else:
                stream.error("%r" % self.tok)
        else:
            stream.take(len(self.tok))
            output.append(self.tok)

    def __str__(self):
        return "%r" % self.tok

    def __repr__(self):
        return "Literal(%r)" % self.tok


class Filler(Parser):
    def __init__(self, tok, no_error=False):
        self.tok = tok
        self.no_error = no_error

    def __call__(self, stream, output):
        if stream.peek(len(self.tok)) != self.tok:
            if self.no_error:
                # Suppress error message.
                # (It's handled by RINGError in error.py)
                stream.error(None)
            else:
                stream.error("%r" % self.tok)
        else:
            stream.take(len(self.tok))

    def __str__(self):
        return "%r" % self.tok

    def __repr__(self):
        return "Literal(%r)" % self.tok


class DeprecatedLiteral(Literal):
    def __call__(self, stream, output):
        import warnings
        Literal.__call__(self, stream, output)
        warnings.warn("Use of '%s' is deprecated in RING inputs" % self.tok)

    def __str__(self):
        return "'%s' (deprecated)" % self.tok

    def __repr__(self):
        return "DeprecatedLiteral(%r)" % self.tok


class Optional(Parser):
    def __init__(self, opt):
        self.opt = opt

    def __call__(self, stream, output):
        with stream:
            stream.parse(self.opt, output)

    def __str__(self):
        return "%s (optional)" % self.opt

    def __repr__(self):
        return "Optional(%r)" % self.opt


class All(Parser):
    def __init__(self, *reqs):
        self.reqs = reqs

    def __call__(self, stream, output):
        for req in self.reqs:
            stream.parse(req, output)

    def __str__(self):
        return '<%s>' % ' '.join('%s' % req for req in self.reqs)

    def __repr__(self):
        return 'All(%s)' % ', '.join('%r' % req for req in self.reqs)


class Either(Parser):
    def __init__(self, *alts):
        self.alts = alts

    def __call__(self, stream, output):
        for alt in self.alts:
            with stream:
                stream.parse(alt, output)
                return
        raise stream.current_error

    def __str__(self):
        return '<%s>' % ' | '.join('%s' % (alt,) for alt in self.alts)

    def __repr__(self):
        return 'Either(%s)' % ', '.join("%r" % (alt,) for alt in self.alts)


class ZeroOrMore(Parser):
    def __init__(self, what):
        self.what = what

    def __call__(self, stream, output):
        with stream:
            stream.parse(self.what, output)
        while not stream.has_error:
            with stream:
                stream.parse(self.what, output)

    def __str__(self):
        return '%s (zero or more)' % self.what

    def __repr__(self):
        return 'ZeroOrMore(%r)' % self.what


class Literals(Either):
    def __init__(self, literal_strings):
        # Sort literals so that parser attempts to match longer literals first.
        # i.e.:  'Cl' should match chlorine and not carbon (with an 'l' left
        # over)
        sorted_strings = sorted(
            literal_strings, key=lambda s: len(s), reverse=True)
        Either.__init__(self, *(Literal(s, no_error=True)
                        for s in sorted_strings))

    def __call__(self, stream, output):
        try:
            return Either.__call__(self, stream, output)
        except RINGSyntaxError:
            if hasattr(self, 'name'):
                stream.error('<' + self.name + '>')
            else:
                raise

# filler is used in parserstate.take() which will move pointers
# beyond the following symbols.


filler = [' ', '\n', '\t']


class ParseState(object):
    """Parse provided stream into an Abstract Syntax Tree."""
    def __init__(self, grammar, stream, debug=False):
        """Parse provided stream into an Abstract Syntax Tree."""
        self.root, self.rules = grammar
        self.stream = stream
        self.debug = debug

        self.stack = []
        self.has_error = False
        self.current_error = None
        self.sidx = 0    # position in the actual string.
        self.lineno = 1  # Human readable pointer for line number
        self.colno = 1   # Human readable pointer for column number
        self.skip_filler()
        if self.debug:
            self.depth = 0

    def error(self, mesg):
        """Raise a RINGSyntaxError exception with the given message.  Line
        and column numbers are provided to the exception object automatically.
        """
        raise RINGSyntaxError(mesg, self.lineno, self.colno, self.stream)

    def peek(self, n=1):
        """Return the next n characters in the character stream.  Do not
        update the stream pointer.
        """
        return self.stream[self.sidx:self.sidx + n]

    def take(self, n=1):
        """Return the next n characters in the character stream and update the
        stream pointer.
        """
        # Move pointer to after the space.
        s = self.stream[self.sidx:self.sidx + n]
        for chr in s:
            if chr == '\n':  # if string contains enter, move to next line.
                self.lineno += 1
                self.colno = 1
            else:
                self.colno += 1
        self.sidx += n
        # skip '\n' and  ' '
        self.skip_filler()

        return s

    def skip_filler(self):
        while self.peek() in filler:
            if self.peek() == '\n':
                self.lineno += 1
                self.colno = 1
                self.sidx += 1
            else:
                self.colno += 1
                self.sidx += 1

    def __enter__(self):
        self.stack.append((self.sidx, self.lineno, self.colno))
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None:
            (self.sidx, self.lineno, self.colno) = self.stack.pop()

            if exc_type is RINGSyntaxError:
                if self.current_error is not None:
                    exc_value.update(self.current_error)
                self.current_error = exc_value
                self.has_error = True
                return True
        else:
            self.stack.pop()
            self.has_error = False
        return False

    def parse(self, what=None, output=None):
        """Parse stream."""
        if what is None:  # Root item.
            output = []
            self.parse(self.root, output)
            return output[0]
        # Call the parse class
        if isinstance(what, Parser):
            inside_output = output[:]
            what(self, inside_output)
            output[:] = inside_output

        # make the what into class, and call the parse on the class
        elif isinstance(what, str):
            inside_output = [RINGToken(what)]
            if self.debug:
                print('%s%r {' % (' '*(self.depth*4), what))
                self.depth += 1
            try:
                self.parse(self.rules[what], inside_output)
            except RINGSyntaxError as exc:
                if self.debug:
                    print('%s%s' % (' '*(self.depth*4), exc))
                raise
            else:
                output.append(inside_output)
                if self.debug:
                    print('%s%r' % (' '*(self.depth*4), output))
            finally:
                if self.debug:
                    self.depth -= 1
                    print('%s} %r' % (' '*(self.depth*4), what))
        else:
            raise TypeError("Don't know how to parse type %s" % (what,))


def parse(stream, strict=False):
    if strict:
        from .Grammar import strict_grammar as grammar
    else:
        from .Grammar import enhanced_grammar as grammar
    return ParseState(grammar, stream).parse()
