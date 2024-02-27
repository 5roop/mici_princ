from pathlib import Path
import pandas as pd
from lxml import etree as ET

template = ET.fromstring(Path("exb_template.xml").read_bytes())
dia = Path("diarization/MP_01.rttm")
ali = Path("kaldi_output/ali.ctm")
diadf = pd.read_csv(
    dia,
    sep="\s",
    names="a b c start duration aa bb speaker aaa bbb".split(),
    usecols="start duration speaker".split(),
)
diadf["end"] = diadf.start + diadf.duration
alidf = pd.read_csv(
    ali,
    delimiter="\s",
    names="sampleid channels start duration word confidence".split(),
)
alidf["end"] = alidf.start + alidf.duration


def assign_speaker_to_alignment(
    alidf: pd.DataFrame, diadf: pd.DataFrame
) -> pd.DataFrame:
    def find_speaker_for_row(row, diadf):
        start = row["start"]
        end = row["end"]
        scores = []
        for i, drow in diadf.iterrows():
            dstart = drow["start"]
            dend = drow["end"]
            overlap = max(0, min(end, dend) - max(start, dstart))
            scores.append(overlap / (end - start))
        diadf["score"] = scores
        speaker = diadf.loc[diadf.score.argmax(), "speaker"]
        del diadf["score"]
        return speaker

    alidf["speaker_name"] = [
        find_speaker_for_row(row, diadf) for i, row in alidf.iterrows()
    ]
    return alidf


newalidf = assign_speaker_to_alignment(alidf, diadf)


def add_df_to_template(exb: ET.Element, df: pd.DataFrame) -> ET.Element:
    """Adds the transcription data from df to the template.

    Right now, only inclusion of one transcription is enabled.

    Args:
        exb (ET.Element): parsed EXB template
        df (pd.DataFrame): dataframe with columns start, end, speaker_name, duration, whisper

    Returns:
        ET.Element: template with transcription tiers.
    """
    df2 = df[["start", "end", "speaker_name", "duration", "word"]].copy()
    df2["speaker_name"] = df.speaker_name
    df2 = df2.rename(columns={"word": "text"})
    df = df2
    df["text"] = df.text.fillna("")
    df["start"] = [float(f"{i:0.3f}") for i in df.start]
    df["end"] = [float(f"{i:0.3f}") for i in df.end]
    # Add speakers:
    for speaker_name in sorted(df.speaker_name.unique()):
        speaker = ET.Element("speaker", attrib={"id": speaker_name})
        abbreviation = ET.Element("abbreviation")
        abbreviation.text = speaker_name
        speaker.append(abbreviation)
        exb.find(".//speakertable").append(speaker)

    # Add <tli>:
    timeline = exb.find(".//common-timeline")
    N = len(timeline.findall(".//tli"))
    for t in sorted(list(set(df.start.values).union(set(df.end.values)))):
        tli = ET.Element("tli", attrib={"id": f"T{N}", "time": str(t)})
        timeline.append(tli)
        N += 1
    # Sort timeline
    timeline[:] = sorted(timeline, key=lambda child: float(child.get("time")))
    # Prepare inverse mapper (seconds -> id):
    mapper = {tli.get("time"): tli.get("id") for tli in timeline.findall("tli")}
    # Add new tier(s):
    for speaker_name in sorted(df.speaker_name.unique()):
        tier = ET.Element(
            "tier",
            attrib=dict(
                id=speaker_name,
                category="v",
                type="t",
                display_name=speaker_name,
                speaker=speaker_name,
            ),
        )
        for i, row in df[df.speaker_name == speaker_name].iterrows():
            event = ET.Element(
                "event",
                attrib=dict(
                    start=mapper.get(str(row["start"])),
                    end=mapper.get(str(row["end"])),
                ),
            )
            event.text = (
                row["text"]
                # if float(row["duration"]) >= min_duration_seconds
                # else "-"
            )
            tier.append(event)
        exb.find("basic-body").append(tier)
    return exb


exb = add_df_to_template(template, newalidf)
exb.find(".//referenced-file").set("url", dia.with_suffix(".wav").name)
ET.indent(exb, space="\t")
out_path = Path("exb_out/") / dia.with_suffix(".exb").name
exb.getroottree().write(
    Path(out_path),
    pretty_print=True,
    encoding="utf8",
    xml_declaration='<?xml version="1.0" encoding="UTF-8"?>',
)

Path(out_path).write_text(
    Path(out_path)
    .read_text()
    .replace(
        "<?xml version='1.0' encoding='UTF8'?>",
        '<?xml version="1.0" encoding="UTF-8"?>',
    )
)
