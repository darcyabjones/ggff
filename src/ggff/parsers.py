#!/usr/bin/env python3

import re

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Union
    from typing import Optional
    from typing import Sequence
    from typing import Mapping, Dict
    from typing import TypeVar
    from typing import TextIO
    from typing import Pattern
    from typing import Callable

    T = TypeVar("T")

from .higher import or_else

MULTISPACE_REGEX = re.compile(r"\s+")


class ValueParseError(Exception):

    def __init__(
        self,
        got: str,
        expected: str,
        template: str = "Could not parse value {got} as {expected}."
    ):
        self.got = got
        self.expected = expected
        self.template = template

    def __str__(self) -> str:
        return self.template.format(got=self.got, expected=self.expected)

    def as_field_error(
        self,
        field: str,
        kind: str = "field"
    ) -> "FieldParseError":
        return FieldParseError(field, str(self), kind=kind)

    def as_line_error(
        self,
    ) -> "LineParseError":
        return LineParseError(str(self))

    def as_block_error(self, line: Optional[int] = None) -> "BlockParseError":
        return BlockParseError(line, str(self))

    def as_parse_error(
        self,
        filename: Optional[str] = None,
        line: Optional[int] = None
    ) -> "ParseError":
        return ParseError(filename, line, str(self))


class FieldParseError(Exception):

    def __init__(self, field: str, message: str, kind: str = "column"):
        self.field = field
        self.message = message
        self.kind = kind
        return

    def __str__(self) -> str:
        return f"In {self.kind} '{self.field}': {self.message}"

    def as_line_error(self) -> "LineParseError":
        return LineParseError(str(self))

    def as_block_error(self, line: Optional[int] = None) -> "BlockParseError":
        return BlockParseError(line, str(self))

    def as_parse_error(
        self,
        filename: Optional[str] = None,
        line: Optional[int] = None,
    ) -> "ParseError":
        return ParseError(filename, line, str(self))


class LineParseError(Exception):
    def __init__(self, message: str):
        self.message = message
        return

    def __str__(self) -> str:
        return self.message

    def as_block_error(self, line: Optional[int] = None) -> "BlockParseError":
        return BlockParseError(line, str(self))

    def as_parse_error(
        self,
        filename: Optional[str] = None,
        line: Optional[int] = None,
    ) -> "ParseError":
        return ParseError(filename, line, str(self))


class BlockParseError(Exception):

    def __init__(self, line: Optional[int], message: str):
        self.line = line
        self.message = message
        return

    def as_block_error(self, line: Optional[int] = None) -> "BlockParseError":
        if (line is None) and (self.line is None):
            nline = None
        else:
            nline = or_else(0, line) + or_else(0, self.line)

        return BlockParseError(
            nline,
            self.message
        )

    def as_parse_error(
        self,
        filename: Optional[str] = None,
        line: Optional[int] = None,
    ) -> "ParseError":
        if (line is None) and (self.line is None):
            line = None
        else:
            line = or_else(0, line) + or_else(0, self.line)

        return ParseError(
            filename,
            line,
            self.message
        )


class ParseError(Exception):
    """ Some aspect of parsing failed. """

    def __init__(
        self,
        filename: Optional[str],
        line: Optional[int],
        message: str
    ):
        self.filename = filename
        self.line = line
        self.message = message
        return

    def add_filename_from_handle(self, handle: TextIO) -> "ParseError":
        if hasattr(handle, "name"):
            self.filename = handle.name
        else:
            self.filename = None

        return self


def split_at_eq(
    fn: Callable[[str], T],
    expected_lhs: str
) -> Callable[[str], Union[T, ValueParseError]]:

    def inner(value: str) -> Union[T, ValueParseError]:
        """ Parse a field of a "key value" type, returning the value. """
        sfield = value.split("=", maxsplit=1)
        if len(sfield) != 2:
            return ValueParseError(
                value,
                f"{expected_lhs}=<value>",
                "Expected a string of the form '{expected}', but got '{got}'",
            )
        elif sfield[0] != expected_lhs:
            return ValueParseError(
                value,
                f"{expected_lhs} <value>",
                "Expected a string of the form '{expected}', but got '{got}'",
            )
        else:
            return fn(sfield[1])

    return inner


def split_at_multispace(
    fn: Callable[[str], T],
    expected_lhs: str
) -> Callable[[str], Union[T, ValueParseError]]:

    def inner(value: str) -> Union[T, ValueParseError]:
        """ Parse a field of a "key value" type, returning the value. """
        sfield = MULTISPACE_REGEX.split(value, maxsplit=1)
        if len(sfield) != 2:
            return ValueParseError(
                value,
                f"{expected_lhs} <value>",
                "Expected a string of the form '{expected}', but got '{got}'",
            )
        elif sfield[0] != expected_lhs:
            return ValueParseError(
                value,
                f"{expected_lhs} <value>",
                "Expected a string of the form '{expected}', but got '{got}'",
            )
        else:
            return fn(sfield[1])

    return inner


def is_value(expected: str) -> Callable[[str], bool]:
    """ """

    def inner(value: str) -> bool:
        return value == expected

    return inner


def is_not_empty_string(value: str) -> str:
    if value.strip() == "":
        raise ValueError("String cannot be empty")
    else:
        return value


def is_one_of(
    options: Sequence[str]
) -> Callable[[str], Union[ValueParseError, str]]:

    def inner(value: str) -> Union[ValueParseError, str]:
        if value in options:
            return value
        else:
            raise ValueParseError(
                value,
                "{" + ", ".join(options) + "}",
                "Invalid value '{got}', must be one of {expected}."
            )

    return inner


def parse_bool(
    true_value: str,
    false_value: str
) -> Callable[[str], Union[ValueParseError, bool]]:
    """ """

    def inner(value: str) -> Union[ValueParseError, bool]:
        if value == true_value:
            return True
        elif value == false_value:
            return False
        else:
            raise ValueParseError(
                value,
                f"{true_value} or {false_value}",
                "Invalid value '{got}', must be either {expected}."
            )

    return inner


def parse_bool_options(
    true_values: Sequence[str],
    false_values: Sequence[str]
) -> Callable[[str], Union[ValueParseError, bool]]:
    """ """
    trues = set(true_values)
    falses = set(false_values)

    def inner(value: str) -> Union[ValueParseError, bool]:
        if value in trues:
            return True
        elif value in falses:
            return False
        else:
            raise ValueParseError(
                value,
                f"{repr(true_values)} or {repr(false_values)}",
                "Invalid value '{got}', must be one of {expected}."
            )

    return inner


def parse_regex(
    regex: Pattern
) -> Callable[[str], Union[ValueParseError, Dict[str, str]]]:

    def inner(value: str) -> Union[ValueParseError, Dict[str, str]]:
        matches = regex.match(value)
        if matches is None:
            return ValueParseError(
                value,
                str(regex),
                ("Value '{got}' does not match the "
                 "regular expression {expected}.")
            )
        else:
            return matches.groupdict()

    return inner


def parse_or_none(
    fn: Callable[[str], Union[ValueParseError, T]],
    none_value: str,
) -> Callable[[str], Union[ValueParseError, T, None]]:
    """ If the value is the same as the none value, will return None.
    Otherwise will attempt to run the fn with field and field name as the
    first and 2nd arguments.
    """

    def inner(value: str) -> Union[ValueParseError, T, None]:
        if value == none_value:
            return None

        val = fn(value)
        if isinstance(val, ValueParseError):
            val.template = (
                f"{val.template}. The value may also be "
                f"'{none_value}', which will be interpreted as None."
            )

        return val

    return inner


def get_from_dict_or_err(
    key: str
) -> Callable[[Mapping[str, T]], Union[ValueParseError, T]]:

    def inner(d: Mapping[str, T]) -> Union[ValueParseError, T]:
        if key in d:
            return d[key]
        else:
            return ValueParseError(key, "", "{got} was missing.")

    return inner


def parse_value(
    fn: Callable[[str], T],
    expected: str,
    etemplate: Optional[str] = None,
) -> Callable[[str], Union[ValueParseError, T]]:

    def inner(value: str) -> Union[ValueParseError, T]:
        try:
            return fn(value)
        except ValueParseError as e:
            return e
        except ValueError:
            if etemplate is not None:
                return ValueParseError(value, expected, etemplate)
            else:
                return ValueParseError(value, expected)

    return inner


parse_int = parse_value(int, "an integer")
parse_float = parse_value(float, "a float")
parse_str = parse_value(is_not_empty_string, "a non-empty string")


def parse_field(
    fn: Callable[[str], Union[ValueParseError, T]],
    field: str,
    kind: str = "field",
) -> Callable[[str], Union[FieldParseError, T]]:

    def inner(value: str) -> Union[FieldParseError, T]:
        result = fn(value)
        if isinstance(result, ValueParseError):
            return result.as_field_error(field, kind=kind)
        else:
            return result

    return inner


def raise_it(fn: Callable[[str], Union[Exception, T]]) -> Callable[[str], T]:

    def inner(value: str) -> T:
        fval = fn(value)
        if isinstance(fval, Exception):
            raise fval
        else:
            return fval

    return inner
