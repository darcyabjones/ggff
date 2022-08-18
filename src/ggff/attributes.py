#!/usr/bin/env python3

import abc
from typing import Mapping, Dict
from typing import Sequence, List
from typing import Iterable, Iterator
from typing import TypeVar
from typing import Union, Optional
from typing import cast

from enum import Enum

from .parsers import (
    raise_it,
    parse_bool_options,
    parse_int,
    parse_field,
    parse_str,
    is_one_of,
    FieldParseError,
)

from .higher import fmap

T = TypeVar('T')

GFF3_KEY_TO_ATTR: Dict[str, str] = {
    "ID": "id",
    "Name": "name",
    "Alias": "alias",
    "Parent": "parent",
    "Target": "target",
    "Gap": "gap",
    "Derives_from": "derives_from",
    "Note": "note",
    "Dbxref": "dbxref",
    "Ontology_term": "ontology_term",
    "Is_circular": "is_circular",
}

GFF3_ATTR_TO_KEY: Dict[str, str] = {v: k for k, v in GFF3_KEY_TO_ATTR.items()}

GFF3_WRITE_ORDER: List[str] = [
    "ID",
    "Name",
    "Alias",
    "Parent",
    "Target",
    "Gap",
    "Derives_from",
    "Note",
    "Dbxref",
    "Ontology_term",
    "Is_circular",
]


attr_is_circular = raise_it(parse_field(
    parse_bool_options(["true", "TRUE", "True"],
                       ["false", "FALSE", "False"]),
    "is_circular",
    "attributes"
))

attr_target_id = raise_it(parse_field(parse_str, "target.id" "attributes"))
attr_target_start = raise_it(parse_field(
    parse_int,
    "target.start",
    "attributes"
))
attr_target_end = raise_it(parse_field(parse_int, "target.end", "attributes"))
attr_target_strand = raise_it(parse_field(
    is_one_of(["+", "-"]),
    "target.strand",
    "attributes"
))
attr_gap_code = raise_it(parse_field(
    is_one_of(["M", "I", "D", "F", "R"]),
    "gap",
    "attributes"
))
attr_gap_len = raise_it(parse_field(parse_int, "gap", "attributes"))


def parse_attr_list(string: str) -> List[str]:
    return list(f.strip() for f in string.strip(", ").split(","))


def flatten_list_of_lists(li: Iterable[Iterable[T]]) -> Iterator[T]:
    for i in li:
        for j in i:
            yield j
    return


def attr_escape(string: str) -> str:
    return (
        string
        .replace("%", "%25")
        .replace(";", "%3B")
        .replace(",", "%2C")
        .replace("=", "%3D")
        .replace("\t", "%09")
        .replace("\n", "%0A")
        .replace("\r", "%0D")
    )


def attr_unescape(string: str) -> str:
    return (
        string
        .replace("%3B", ";")
        .replace("%2C", ",")
        .replace("%3D", "=")
        .replace("%09", "\t")
        .replace("%0A", "\n")
        .replace("%0D", "\r")
        .replace("%25", "%")
    )


class TargetStrand(Enum):
    PLUS = 0
    MINUS = 1

    def __str__(self) -> str:
        into_str_map = ["+", "-"]
        return into_str_map[self.value]

    def __repr__(self) -> str:
        return "TargetStrand.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "TargetStrand":
        from_str_map = {
            "+": cls.PLUS,
            "-": cls.MINUS,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")


class Target(object):

    """
    Indicates the target of a nucleotide-to-nucleotide or
    protein-to-nucleotide alignment.
    The format of the value is "target_id start end [strand]",
    where strand is optional and may be "+" or "-".
    If the target_id contains spaces, they must be escaped as hex escape %20.
    """

    def __init__(
        self,
        target_id: str,
        start: int,
        end: int,
        strand: Optional[TargetStrand] = None,
    ) -> None:
        self.target_id = target_id
        self.start = start
        self.end = end
        self.strand = strand
        return

    def __repr__(self) -> str:
        if self.strand is None:
            return f"Target('{self.target_id}', {self.start}, {self.end})"
        else:
            strand = repr(self.strand)
            return (
                f"Target('{self.target_id}', {self.start}, "
                f"{self.end}, {strand})"
            )

    def __str__(self) -> str:
        target_id = attr_escape(self.target_id).replace(" ", "%20")

        # Recode back to 1 based (inclusive)
        start = self.start + 1

        if self.strand is None:
            return "{} {} {}".format(target_id, start, self.end)
        else:
            return "{} {} {} {}".format(
                target_id,
                start,
                self.end,
                self.strand
            )

    @classmethod
    def parse(cls, string: str, unescape: bool = True) -> "Target":
        split_string = string.strip().split(" ")

        if len(split_string) < 3:
            raise FieldParseError(
                "target",
                ("Too few fields in Target string. Expected 3 or 4, "
                 f"but got {len(split_string)}."),
                "attribute"
            )
        elif len(split_string) > 4:
            raise FieldParseError(
                "target",
                ("Too many fields in Target strand. Expected 3 or 4, "
                 f"but got {len(split_string)}."),
                "attribute"
            )

        if unescape:
            target_id = attr_target_id(
                attr_unescape(split_string[0])
                .replace("%20", " ")
            )
        else:
            target_id = attr_target_id(split_string[0])
        # We want 0 based exclusive
        start = attr_target_start(split_string[1]) - 1
        end = attr_target_end(split_string[2])

        if len(split_string) == 3:
            return cls(target_id, start, end)
        else:
            strand = TargetStrand.parse(attr_target_strand(split_string[3]))
            return cls(target_id, start, end, strand)


class GapCode(Enum):

    MATCH = 0
    INSERT = 1
    DELETE = 2
    FORWARD = 3
    REVERSE = 4

    def __str__(self) -> str:
        into_str_map = ["M", "I", "D", "F", "R"]
        return into_str_map[self.value]

    def __repr__(self) -> str:
        return f"GapCode.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "GapCode":
        from_str_map = {
            "M": cls.MATCH,
            "I": cls.INSERT,
            "D": cls.DELETE,
            "F": cls.FORWARD,
            "R": cls.REVERSE,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")
        return


class GapElement(object):

    def __init__(self, code: GapCode, length: int) -> None:
        """ """
        assert length >= 0

        self.code = code
        self.length = length
        return

    def __str__(self) -> str:
        return "{}{}".format(self.code, self.length)

    def __repr__(self) -> str:
        code = repr(self.code)
        return f"GapElement({code}, {self.length})"

    @classmethod
    def parse(cls, string: str) -> "GapElement":
        string = string.strip()
        code = GapCode.parse(attr_gap_code(string[:1]))
        length = attr_gap_len(string[1:])
        return cls(code, length)


class Gap(object):

    def __init__(self, elements: Sequence[GapElement]) -> None:
        self.elements = list(elements)
        return

    def __str__(self) -> str:
        return " ".join(map(str, self.elements))

    def __repr__(self) -> str:
        elems = repr(list(self.elements))
        return f"Gap({elems})"

    @classmethod
    def parse(cls, string: str) -> "Gap":
        split_string = string.strip().split(" ")
        elements = [GapElement.parse(s) for s in split_string if s != '']
        return cls(elements)


# This is necessary until Self becomes standard
AttrType = TypeVar("AttrType", bound="Attributes")


class Attributes(abc.ABC):

    @abc.abstractmethod
    def as_str(self, escape: bool = True) -> str:
        raise NotImplementedError()

    def __str__(self) -> str:
        return self.as_str()

    @abc.abstractmethod
    def __repr__(self) -> str:
        ...

    @classmethod
    @abc.abstractmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = True,
    ) -> AttrType:
        ...

    @property
    @abc.abstractmethod
    def id(self) -> str:
        ...

    @id.setter
    @abc.abstractmethod
    def id(self, val: str):
        ...

    @property
    @abc.abstractmethod
    def parent(self) -> List[str]:
        ...

    @parent.setter
    @abc.abstractmethod
    def parent(self, vals: Sequence[str]):
        ...


GFF3AttributesType = TypeVar("GFF3AttributesType", bound="GFF3Attributes")


class GFF3Attributes(Attributes):

    def __init__(
        self,
        id: Optional[str] = None,
        name: Optional[str] = None,
        alias: Optional[Sequence[str]] = None,
        parent: Optional[Sequence[str]] = None,
        target: Optional[Target] = None,
        gap: Optional[Gap] = None,
        derives_from: Optional[Sequence[str]] = None,
        note: Optional[Sequence[str]] = None,
        dbxref: Optional[Sequence[str]] = None,
        ontology_term: Optional[Sequence[str]] = None,
        is_circular: bool = False,
        custom: Optional[Mapping[str, Sequence[str]]] = None,
    ) -> None:
        self.id = id
        self.name = name

        if alias is not None:
            self.alias: List[str] = list(alias)
        else:
            self.alias = []

        if parent is not None:
            self.parent: List[str] = list(parent)
        else:
            self.parent = []

        self.target = target
        self.gap = gap

        if derives_from is not None:
            self.derives_from: List[str] = list(derives_from)
        else:
            self.derives_from = []

        if note is not None:
            self.note: List[str] = list(note)
        else:
            self.note = []

        if dbxref is not None:
            self.dbxref: List[str] = list(dbxref)
        else:
            self.dbxref = []

        if ontology_term is not None:
            self.ontology_term: List[str] = list(ontology_term)
        else:
            self.ontology_term = []

        self.is_circular = is_circular

        self.custom: Dict[str, Union[str, List[str]]] = {}
        if custom is not None:
            for k, v in custom.items():
                if isinstance(v, str):
                    self.custom[k] = v
                else:
                    self.custom[k] = list(v)

        return

    @classmethod
    def _parse_list_of_strings(
        cls,
        value: str,
        strip_quote: bool = False,
        unescape: bool = True
    ) -> List[str]:
        """ Parses a gff attribute list of strings. """
        if value == "":
            return []

        if strip_quote:
            value = value.strip("\"' ")

        if unescape:
            return [attr_unescape(v) for v in parse_attr_list(value)]
        else:
            return parse_attr_list(value)

    @classmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = True,
    ) -> GFF3AttributesType:
        if string.strip() in (".", ""):
            return cls()

        fields = (
            f.split("=", maxsplit=1)
            for f
            in string.strip(" ;").split(";")
        )

        if strip_quote:
            kvpairs: Dict[str, str] = {
                k.strip(): v.strip(" '\"")
                for k, v
                in fields
            }
        else:
            kvpairs = {
                k.strip(): v.strip()
                for k, v
                in fields
            }

        if unescape:
            id = fmap(attr_unescape, kvpairs.pop("ID", None))
        else:
            id = kvpairs.pop("ID", None)

        if id == "":
            id = None

        if unescape:
            name = fmap(attr_unescape, kvpairs.pop("Name", None))
        else:
            name = kvpairs.pop("Name", None)

        if name == "":
            name = None

        alias = cls._parse_list_of_strings(
            kvpairs.pop("Alias", ""),
            strip_quote,
            unescape
        )

        parent = cls._parse_list_of_strings(
            kvpairs.pop("Parent", ""),
            strip_quote,
            unescape
        )

        target: Optional[Target] = fmap(
            lambda x: Target.parse(x, unescape),
            kvpairs.pop("Target", None)
        )

        gap = fmap(Gap.parse, kvpairs.pop("Gap", None))

        derives_from = cls._parse_list_of_strings(
            kvpairs.pop("Derives_from", ""),
            strip_quote,
            unescape
        )

        note = cls._parse_list_of_strings(
            kvpairs.pop("Note", ""),
            strip_quote,
            unescape
        )

        dbxref = cls._parse_list_of_strings(
            kvpairs.pop("Dbxref", ""),
            strip_quote,
            unescape
        )

        ontology_term = cls._parse_list_of_strings(
            kvpairs.pop("Ontology_term", ""),
            strip_quote,
            unescape
        )

        is_circular = attr_is_circular(kvpairs.pop("Is_circular", "false"))

        custom: dict[str, Union[str, List[str]]] = dict()
        for k, v in kvpairs.items():
            if "," in v:
                custom[k] = cls._parse_list_of_strings(
                    v, strip_quote, unescape)
            elif v != "":
                custom[k] = v

        return cls(
            id,
            name,
            alias,
            parent,
            target,
            gap,
            derives_from,
            note,
            dbxref,
            ontology_term,
            is_circular,
            custom
        )

    def is_empty(self) -> bool:  # noqa
        # Yes, this could be written as single boolean comparison.
        # But it gets so long that it's hard to understand.
        if len(self.custom) > 0:
            return False
        elif self.id is not None:
            return False
        elif self.name is not None:
            return False
        elif len(self.alias) > 0:
            return False
        elif len(self.parent) > 0:
            return False
        elif self.target is not None:
            return False
        elif self.gap is not None:
            return False
        elif len(self.derives_from) > 0:
            return False
        elif len(self.note) > 0:
            return False
        elif len(self.dbxref) > 0:
            return False
        elif len(self.ontology_term) > 0:
            return False
        elif self.is_circular:
            return False
        else:
            return True

    def __str__(self) -> str:
        return self.as_str()

    def as_str(self, escape: bool = True) -> str:
        # Avoid having an empty string at the end.
        if self.is_empty():
            return "."

        keys = []
        keys.extend(GFF3_WRITE_ORDER)
        keys.extend(self.custom.keys())

        kvpairs = []
        for key in keys:
            value = self[key]
            if value is None or value == []:
                continue
            elif key == "Is_circular" and not value:
                continue

            if escape:
                key = attr_escape(key)

            if isinstance(value, list) and escape:
                value = ",".join(attr_escape(str(v)) for v in value)

            elif isinstance(value, list):
                value = ",".join(str(v) for v in value)

            elif isinstance(value, bool):
                value = "true" if value else "false"

            else:
                if escape:
                    value = attr_escape(str(value))
                else:
                    value = str(value)

            kvpairs.append((key, value))

        return ";".join(f"{k}={v}" for k, v in kvpairs)

    def __repr__(self) -> str:
        param_names = [GFF3_KEY_TO_ATTR[k] for k in GFF3_WRITE_ORDER]
        param_names.append("custom")

        parameters = []
        for param in param_names:
            value = getattr(self, param)

            if value is None or value == []:
                continue

            if isinstance(value, list):
                value = repr(list(value))
            else:
                value = repr(value)

            parameters.append(f"{param}={value}")

        joined_parameters = ", ".join(parameters)
        return f"GFF3Attributes({joined_parameters})"

    def __getitem__(
        self,
        key: str,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        if key in GFF3_KEY_TO_ATTR:
            return getattr(self, GFF3_KEY_TO_ATTR[key])
        else:
            return self.custom[key]

    def __setitem__(
        self,
        key: str,
        value: Union[str, Sequence[str], Target, Gap, bool, None],
    ) -> None:
        """ If the key is an attr set it, otherwise update the custom dict."""

        if key in GFF3_KEY_TO_ATTR:
            setattr(self, GFF3_KEY_TO_ATTR[key], value)
        else:
            self.custom[key] = cast(str, value)  # Cast is for mypy
        return

    def __delitem__(self, key: str) -> None:
        """ Removes an item from the custom dictionary or resets
        attribute to default """

        if key in ("ID", "Name", "Target", "Gap"):
            setattr(self, GFF3_KEY_TO_ATTR[key], None)

        elif key in ("Alias", "Parent", "Derives_from", "Note",
                     "Dbxref", "Ontology_term"):
            setattr(self, GFF3_KEY_TO_ATTR[key], [])

        elif key == "Is_circular":
            setattr(self, GFF3_KEY_TO_ATTR[key], False)

        else:
            del self.custom[key]

        return

    def get(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        """ Gets an atrribute or element from the custom dict. """

        if key in GFF3_KEY_TO_ATTR:
            return getattr(self, GFF3_KEY_TO_ATTR[key])
        else:
            return self.custom.get(key, default)

    def pop(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        """ Removes an item from the attributes and returns its value.

        If the item is one of the reserved GFF3 terms, the
        value is reset to the default.
        """

        if key in GFF3_KEY_TO_ATTR:
            value = getattr(self, GFF3_KEY_TO_ATTR[key])
            del self[GFF3_KEY_TO_ATTR[key]]
            return value

        else:
            return self.custom.pop(key, default)
