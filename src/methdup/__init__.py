"""Public API for the methdup BAM deduplicator."""

from .cli import main, run
from .models import Counters, DeduplicationError
from .processor import deduplicate

__all__ = ["Counters", "DeduplicationError", "deduplicate", "main", "run"]
