from pathlib import Path
import pandas as pd

ali = Path("output/ali.ctm")
df = pd.read_csv(ali, delimiter="\s", names="sampleid channels start duration word confidence".split())
df.to_csv("002_aligned.csv", index=False)