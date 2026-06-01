"""Backward-compatible entry point for the Tahoe-100M embedding pipeline.

New code lives in src/tahoe_pipeline.py. This wrapper keeps the original
filename usable for demos or interview walkthroughs.
"""

from src.tahoe_pipeline import main


if __name__ == "__main__":
    main()
