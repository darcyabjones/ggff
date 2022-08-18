#!/usr/bin/env python3

from typing import TypeVar, Generic
from typing import Sequence, List
from typing import Iterable
from typing import Optional

from .record import GFFRecord

T = TypeVar('T')


class GFF(Generic[T]):

    def __init__(self, records: Optional[Sequence[GFFRecord[T]]] = None):
        self.ids: dict[GFFRecord[T], str] = {}

        if records is None:
            self.records: List[GFFRecord[T]] = []
        else:
            self.records: List[GFFRecord[T]] = list(records)
        return

    def add_record(self, record: GFFRecord[T], force: bool = False):
        return

    def add_records(
        self,
        records: Iterable[GFFRecord[T]],
        force: bool = False
    ):
        return
