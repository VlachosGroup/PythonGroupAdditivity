# -*- coding: utf-8 -*-
from .. Error import RINGError


class Reader(object):
    """
    Reader reads the parsed RING input, and returns the RDkit wrapper objects
    in pgradd.RDkitWrapper.

    Attributes
    ----------
    ast : abstract syntax tree obtrained from parser
    """
    def __init__(self, ast):
        # ast = Abstract Syntax Tree
        self.ast = ast

    def ReadRINGInput(self, tree):
        # Check the type of input
        assert tree[0][0].name in ('Fragment',
                                   'ReactionRule',
                                   'EnumerationQuery')
        # if fragment, molquery is returned
        if tree[0][0].name == 'Fragment':
            from . MolQueryRead import MolQueryReader
            self.type = 'MolQuery'
            return MolQueryReader(tree[0][1:]).Read()
        # if reaction rule, reacitonquery is returned
        elif tree[0][0].name == 'ReactionRule':
            from . ReactionQueryRead import ReactionQueryReader
            self.type = 'ReactionQuery'
            return ReactionQueryReader(tree[0][1:]).Read()
        # TODO enumeration query
        elif tree[0][0].name == 'EnumerationQuery':
            raise NotImplementedError('Coming soon')

    def Read(self):
        # Root tree reading. Check if the input is RINGinput
        assert self.ast[0].name == "RINGInput"
        return self.ReadRINGInput(self.ast[1:])


def Read(text, strict=False):
    """
    Return MolQuery, ReactionQuery, or ReactionNetworkEnumerationQuery by
    interpretting RING input string.

    Parameters
    ----------
    text : string
        Specify string describing chemical structure, elementary reaction, or
        reaction network enumeration rules in RING notation.
    strict : boolean, optional
        If True, then disable use of syntactic extensions such as support for
        "radical electrons".

    Returns
    -------
    Returns RDkit wrapped queries that extends RDkit's functionality:
    * MolQuery if fragment is given in string
    * ReactionQuery if reaction rule is given in string
    * ReactionNetworkEnumerationQuery if enumeration query is given in string

    Raises
    ------
    msr.error.RINGSyntaxError
        If `text` does not conform to RING syntax.
    msr.error.RINGReaderError
        If `text` is invalid RING for non-syntactic reasons.
    """
    from . import Parser
    try:
        return Reader(Parser.parse(text)).Read()
    except RINGError as exc:
        raise exc
