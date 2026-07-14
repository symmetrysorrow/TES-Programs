"""Python interface for the dependency-free PHITS dumpall.dat converter."""

from .converter import ConversionResult, convert

__all__ = ["ConversionResult", "convert"]
