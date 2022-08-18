#!/usr/bin/env python3

from typing import Optional, Union
from typing import Set, List, Dict
from typing import Sequence, Mapping
from typing import Iterable, Iterator
from typing import TypeVar, Generic, Type
from typing import cast
from typing import TextIO

from enum import Enum
from collections import deque

from dataclasses import field

from .higher import fmap
from .parsers import (
    ParseError,
    LineParseError,
    FieldParseError,
    raise_it,
    parse_field,
    parse_int,
    parse_float,
    parse_or_none,
    is_one_of,
    parse_str,
    parse_bool_options,
)
from .attributes import GFF3Attributes
from .gff import GFF

T = TypeVar('T')

rec_seqid = raise_it(parse_field(parse_str, "seqid"))
rec_source = raise_it(parse_field(parse_str, "source"))
rec_type = raise_it(parse_field(parse_str, "type"))
rec_start = raise_it(parse_field(parse_int, "start"))
rec_end = raise_it(parse_field(parse_int, "end"))
rec_score = raise_it(parse_field(parse_or_none(parse_float, "."), "score"))
rec_strand = raise_it(parse_field(is_one_of(["-", "+", ".", "?"]), "strand"))
rec_phase = raise_it(parse_field(is_one_of(["0", "1", "2", "."]), "phase"))


class Strand(Enum):
    PLUS = 0
    MINUS = 1
    UNSTRANDED = 2
    UNKNOWN = 3

    def __str__(self):
        into_str_map: List[str] = ["+", "-", ".", "?"]
        return into_str_map[self.value]

    def __repr__(self):
        return f"Strand.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "Strand":
        from_str_map: Dict[str, Strand] = {
            "+": cls.PLUS,
            "-": cls.MINUS,
            ".": cls.UNSTRANDED,
            "?": cls.UNKNOWN,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")


class Phase(Enum):
    FIRST = 0
    SECOND = 1
    THIRD = 2
    NOT_CDS = 3

    def __str__(self):
        into_str_map: List[str] = ["0", "1", "2", "."]
        return into_str_map[self.value]

    def __repr__(self):
        return f"Phase.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "Phase":
        from_str_map: Dict[str, Phase] = {
            "0": cls.FIRST,
            "1": cls.SECOND,
            "2": cls.THIRD,
            ".": cls.NOT_CDS,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")


class Record(Generic[T]):

    columns: List[str] = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]

    __attribute_cls: Type[T]

    def __init__(
        self,
        seqid: str,
        source: str,
        type: str,
        start: int,
        end: int,
        score: Optional[float] = None,
        strand: Strand = Strand.UNSTRANDED,
        phase: Phase = Phase.NOT_CDS,
        attributes: Optional[T] = None,
        gff: Optional[GFF] = None
    ) -> None:
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase

        if attributes is None:
            self.attributes = self.__attribute_cls()
        else:
            self.attributes = attributes

        self.gff = gff
        return

    def __str__(self) -> str:
        return self.as_str()

    def as_str(self, escape: bool = True) -> str:
        values = []
        for name in self.columns:
            value = getattr(self, name)

            if value is None:
                values.append(".")
            # Convert back to 1-based inclusive indexing.
            elif name == "start":
                values.append(str(value + 1))
            elif name == "attributes" and not escape:
                values.append(value.as_str(escape=False))
            else:
                values.append(str(value))

        return "\t".join(values)

    def __repr__(self) -> str:
        parameters = []
        for col in self.columns:
            val = repr(getattr(self, col))
            parameters.append(f"{val}")

        joined_parameters = ", ".join(parameters)
        return f"GFFRecord({joined_parameters})"

    def length(self) -> int:
        """ Returns the distance between start and end. """
        return self.end - self.start

    def __len__(self) -> int:
        return self.length()

    def add_child(self, child: "GFFRecord") -> None:
        if child not in self.children:
            self.children.append(child)
        if self not in child.parents:
            child.parents.append(self)
        return

    def add_parent(self, parent: "GFFRecord") -> None:
        if parent not in self.parents:
            self.parents.append(parent)

        if self not in parent.children:
            parent.children.append(self)
        return

    def add_children(self, children: Sequence["GFFRecord"]) -> None:
        for child in children:
            self.add_child(child)
        return

    def add_parents(self, parents: Sequence["GFFRecord"]) -> None:
        for parent in parents:
            self.add_parent(parent)
        return

    def add_derivative(self, derivative: "GFFRecord") -> None:
        if derivative not in self.derivatives:
            self.derivatives.append(derivative)
        if self not in derivative.parents:
            derivative.derives_from.append(self)
        return

    def add_derives_from(self, derives_from: "GFFRecord") -> None:
        if derives_from not in self.derives_from:
            self.derives_from.append(derives_from)

        if self not in derives_from.derivatives:
            derives_from.derivatives.append(self)
        return

    def add_derivatives(self, derivatives: Sequence["GFFRecord"]) -> None:
        for derivative in derivatives:
            self.add_derivative(derivative)
        return

    def add_derives_froms(self, derives_froms: Sequence["GFFRecord"]) -> None:
        for derives_from in derives_froms:
            self.add_derives_from(derives_from)
        return

    def update_parents(self) -> None:
        parent_ids = []
        for parent in self.parents:
            parent_id = parent.attributes.id
            assert parent_id is not None
            parent_ids.append(parent_id)

        self.attributes.parent = parent_ids
        return

    def update_derives_from(self) -> None:
        df_ids = []
        for df in self.derives_from:
            df_id = df.attributes.id
            assert df_id is not None
            df_ids.append(df_id)

        self.attributes.derives_from = df_ids
        return

    def traverse_children(
        self,
        sort: bool = False,
        breadth: bool = False,
    ) -> Iterator["GFFRecord"]:
        """ A graph traversal of this `GFFRecord`s children.

        Keyword arguments:
        sort -- Sort the children so that they appear in increasing order.
        breadth -- Perform a breadth first search rather than the default
            depth first search.

        Returns:
        A generator yielding `GFFRecord`s.
        """

        should_reverse = not breadth

        seen: Set[GFFRecord] = set()

        to_visit = deque([self])

        while len(to_visit) > 0:

            if breadth:
                node = to_visit.popleft()
            else:
                node = to_visit.pop()

            # NB uses id() for hashing
            if node in seen:
                continue
            else:
                seen.add(node)

            children = list(node.children)
            if sort:
                children.sort(
                    key=lambda f: (f.seqid, f.start, f.end, f.type),
                    reverse=should_reverse
                )

            to_visit.extend(children)
            yield node

        return None

    def traverse_parents(
        self,
        sort: bool = False,
        breadth: bool = False,
    ) -> Iterator["GFFRecord"]:
        """ A graph traversal of this `GFFRecord`s parents.

        Keyword arguments:
        sort -- Sort the parents so that they appear in increasing order.
        breadth -- Perform a breadth first search rather than the default
            depth first search.

        Returns:
        A generator yielding `GFFRecord`s.
        """

        should_reverse = not breadth

        seen: Set[GFFRecord] = set()
        to_visit = deque([self])

        while len(to_visit) > 0:

            if breadth:
                node = to_visit.popleft()
            else:
                node = to_visit.pop()

            if node in seen:
                continue
            else:
                seen.add(node)

            parents = list(node.parents)
            if sort:
                parents.sort(
                    key=lambda f: (f.seqid, f.start, f.end, f.type),
                    reverse=should_reverse
                )

            to_visit.extend(parents)
            yield node

        return None

    @classmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = True,
    ) -> "GFFRecord":
        """ Parse a gff line string as a `GFFRecord`.

        Keyword arguments:
        string -- The gff line to parse.
        format -- What format the gff file is in.
            Currently only GFF3 is supported.
        strip_quote -- Strip quotes from attributes values. The specification
            says that they should not be stripped, so we don't by default.
        unescape -- Unescape reserved characters in the attributes to their
            original values. I.E. some commas, semicolons, newlines etc.

        Returns:
        A `GFFRecord`
        """

        sline = string.strip().split("\t")
        sline_len = len(sline)
        columns_len = len(cls.columns)

        if sline_len < (columns_len - 1):
            raise ValueError((
                "Line has too few columns. "
                f"Expected: {columns_len}, Encountered: {sline_len}"
            ))

        fields: Dict[str, str] = dict(zip(cls.columns, sline))
        if sline_len == columns_len - 1:
            fields["attributes"] = ""

        # 0-based indexing exclusive
        start = rec_start(fields["start"]) - 1
        end = rec_end(fields["end"])

        if start > end:
            tmp = start
            start = end
            end = tmp
            del tmp

        score = rec_score(fields["score"])
        strand = Strand.parse(rec_strand(fields["strand"]))
        phase = Phase.parse(rec_phase(fields["phase"]))

        attributes = cls.__attribute_cls.parse(
            fields["attributes"],
            strip_quote=strip_quote,
            unescape=unescape,
        )

        return cls(
            rec_seqid(fields["seqid"]),
            rec_source(fields["source"]),
            rec_type(fields["type"]),
            start,
            end,
            score,
            strand,
            phase,
            attributes
        )

    def trim_ends(self, length: int) -> None:
        """ Trim a specified number of bases from each end of the feature. """

        from math import ceil

        if self.length() <= 2:
            length = 0
        elif self.length() < (2 * length):
            length = ceil(self.length() / 4)

        self.start += length
        self.end -= length
        return

    def pad_ends(self, length: int, max_value: Optional[int] = None) -> None:

        if (self.start - length) < 0:
            self.start = 0
        else:
            self.start -= length

        if max_value is not None and (self.end + length) > max_value:
            self.end = max_value
        else:
            self.end += length
        return

    def expand_to_children(self) -> None:
        if len(self.children) == 0:
            return

        min_ = min(c.start for c in self.children)
        max_ = max(c.end for c in self.children)

        if min_ < self.start:
            self.start = min_

        if max_ > self.end:
            self.end = max_
        return

    def shrink_to_children(self) -> None:
        if len(self.children) == 0:
            return

        min_ = min(c.start for c in self.children)
        max_ = max(c.end for c in self.children)

        if min_ > self.start:
            self.start = min_

        if max_ < self.end:
            self.end = max_
        return

    def copy(self) -> "GFFRecord":
        """ You'll still need to update the ID """
        from copy import copy, deepcopy

        node_copy = copy(self)
        node_copy.attributes = deepcopy(self.attributes)

        if node_copy.attributes is not None:
            node_copy.attributes.parent = []

        node_copy.children = []
        node_copy.parents = []
        return node_copy

    @classmethod  # noqa
    def from_file(
        cls,
        handle: TextIO,
        strip_quote: bool = False,
        unescape: bool = True,
    ) -> Iterator["GFFRecord"]:
        """
        Yes yes I need to make this more modular... not time right now.
        """

        from collections import defaultdict

        id_to_record = defaultdict(list)

        to_yield = list()
        lonely_children = set()
        lonely_derivatives = set()

        # Children that are encountered before their parents
        undefined_parents: Dict[str, List[GFFRecord]] = defaultdict(list)

        # derivatives encountered before their parents.
        undefined_derives_froms: Dict[str, List[GFFRecord]] = defaultdict(list)

        for i, line in enumerate(handle):
            if line.startswith("#"):
                continue
            elif line.strip() == "":
                continue

            try:
                record = GFFRecord.parse(line, strip_quote, unescape)
            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)

            id_ = record.attributes.id

            if id_ is not None:
                id_to_record[id_].append(record)

                uparents = undefined_parents.pop(id_, [])
                uderives = undefined_derives_froms.pop(id_, [])
                record.add_children(uparents)
                record.add_derivatives(uderives)

                for r in uparents:
                    lonely_children.discard(r)
                    if r not in lonely_derivatives:
                        to_yield.append(r)

                for r in uderives:
                    lonely_derivatives.discard(r)
                    if r not in lonely_children:
                        to_yield.append(r)

            is_missing_parent = False
            for parent in record.attributes.parent:
                if parent not in id_to_record:
                    undefined_parents[parent].append(record)
                    lonely_children.add(record)
                    is_missing_parent = True
                else:
                    record.add_parents(id_to_record.get(parent, []))

            is_missing_derives_from = False
            for derives_from in record.attributes.derives_from:
                if derives_from not in id_to_record:
                    undefined_derives_froms[derives_from].append(record)
                    lonely_derivatives.add(record)
                    is_missing_derives_from = True
                else:
                    record.add_derives_froms(
                        id_to_record.get(derives_from, []))

            if not (is_missing_parent or is_missing_derives_from):
                to_yield.append(record)

            while len(to_yield) > 0:
                yield to_yield.pop()

        if (len(undefined_parents) > 0) or (len(undefined_derives_froms) > 0):
            upar = set(undefined_parents.keys())
            uder = set(undefined_derives_froms.keys())

            message = [
                "Reached the end of GFF file with undefined references."
            ]

            if len(upar) > 0:
                message.append(f"Expected to find parent: {repr(upar)}.")

            if len(uder) > 0:
                message.append(f"Expected to find derives_from: {repr(uder)}.")

            raise (
                ParseError(filename=None, line=i, message=" ".join(message))
                .add_filename_from_handle(handle)
            )
        return
