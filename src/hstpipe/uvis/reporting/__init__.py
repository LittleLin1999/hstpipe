"""Self-contained UVIS HTML reports; science stages do not depend on rendering."""
from .report import generate_report

__all__ = ["generate_report"]
