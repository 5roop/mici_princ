try:
    exb = snakemake.input.exb
    transcript = snakemake.input.transcript
    jsonfile = snakemake.input.jsonfile
    gold = snakemake.input.gold
    outfile = snakemake.output.out
except Exception as e:
    print(e)
    exb = "/home/rupnik/mici_princ/inspected_exbs_fixed/MP_03.exb"
    transcript = "/home/nikola/transfer/mici_princ/MP_03.txt.txt"
    gold = "/home/nikola/transfer/mici_princ/MP_03.txt"
    jsonfile = "/home/nikola/transfer/mici_princ/MP_03.txt.json"
    outfile = "aligned_jsons/MP_03.json"
from pathlib import Path
from lxml import etree as ET
import json

gold = Path(gold).read_text()
exb = ET.fromstring(Path(exb).read_bytes())
goldenjson = json.loads(Path(jsonfile).read_text())
timeline = {i.get("id"): i.get("time") for i in exb.findall(".//tli")}
events = exb.findall(".//event")
for event in events:
    event.set("start_s", timeline[event.get("start")])
    event.set("end_s", timeline[event.get("end")])

exbevents = sorted(
    [
        [
            i.text.strip(),
            i.get("start_s"),
            i.get("end_s"),
            i.getparent().get("display-name"),
        ]
        for i in events
        if i.text is not None
    ],
    key=lambda l: float(l[1]),
)
assert len(goldenjson) == len(exbevents)
if Path(jsonfile).name == "MP_00.txt.json":
    for i, entry in enumerate(goldenjson):
        word, char_s, char_e = entry
        if char_s > 383:
            char_s -= 1
            char_e -= 1
        if char_e == 383:
            char_e = 382

        goldenjson[i] = [word, char_s, char_e]

final = [i + j for i, j in zip(goldenjson, exbevents)]
import pandas as pd

df = pd.DataFrame(
    final, columns="gold_text char_s char_e text time_s time_e speaker".split()
)
df["time_s_float"] = df.time_s.astype(float)
df = df.sort_values(by="time_s_float")
df["speaker_changed_now"] = df.speaker != df.speaker.shift(1)
df["speaker_will_change"] = df.speaker != df.speaker.shift(-1)
df["time_s"] = df.time_s.astype(float)
df["time_e"] = df.time_e.astype(float)

utterances = list()
current = dict()
utterance = ""
for i, row in df.iterrows():
    is_last = row["speaker_will_change"]
    is_first = row["speaker_changed_now"]
    if is_first & is_last:
        utterance = ""
        current = {k: row[k] for k in "char_s char_e time_s time_e speaker".split()}
        current["words"] = [{k: row[k] for k in "char_s char_e time_s time_e".split()}]
        pad = 0
        while True:
            try:
                if gold[row["char_e"] + pad] in "\n ":
                    break
                else:
                    pad += 1
            except IndexError:
                break
        current["text"] = gold[current["char_s"] : current["char_e"] + pad]
        utterances.append(current)
        print(Fore.GREEN + row["text"])
        print(Fore.BLUE + current["text"])
    elif is_first:
        utterance = ""
        current = dict()
        current["words"] = []
        current["char_s"] = row["char_s"]
        current["time_s"] = row["time_s"]
        current["words"].append(
            {
                "time_s": row["time_s"],
                "time_e": row["time_e"],
                "char_s": row["char_s"],
                "char_e": row["char_e"],
            }
        )
        utterance = utterance + " " + row["text"]
    elif is_last:
        pad = 0
        while True:
            try:
                if gold[row["char_e"] + pad] in "\n ":
                    break
                else:
                    pad += 1
            except IndexError:
                break
        current["char_e"] = row["char_e"]
        current["text"] = gold[current["char_s"] : current["char_e"] + pad]

        current["words"].append(
            {šipak
                "time_s": row["time_s"],
                "time_e": row["time_e"],
                "char_s": row["char_s"],
                "char_e": row["char_e"],  # + pad,
            }
        )
        current["speaker"] = row["speaker"]
        utterances.append(current)
        utterance = utterance + " " + row["text"]
        from colorama import Fore, Back, Style

        print(Fore.GREEN + utterance)
        print(Fore.BLUE + current["text"])
    else:
        utterance = utterance + " " + row["text"]
        current["words"].append(
            {
                "time_s": row["time_s"],
                "time_e": row["time_e"],
                "char_s": row["char_s"],
                "char_e": row["char_e"],
            }
        )
df["original_text"] = df.apply(lambda row: gold[row["char_s"] : row["char_e"]], axis=1)
Path(outfile).write_text(
    json.dumps(utterances, indent=4, ensure_ascii=False, sort_keys=True)
)
