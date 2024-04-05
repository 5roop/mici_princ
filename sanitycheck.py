try:
    exb = snakemake.input.exb
    transcript = snakemake.input.transcript
    jsonfile = snakemake.input.jsonfile
    outfile = snakemake.output.out
except:
    exb = "/home/nikola/transfer/mici_princ/MP_25.exb"
    transcript = "/home/nikola/transfer/mici_princ/MP_25.txt.txt"
    jsonfile = "/home/nikola/transfer/mici_princ/MP_25.txt.json"
    outfile = "inspected_exbs_fixed/MP_25.exb"
from pathlib import Path
from lxml import etree as ET

exb = ET.fromstring(Path(exb).read_bytes())

for tier in exb.findall(".//tier"):
    tier.set(
        "display-name",
        tier.get("display-name")
        .replace(" [v]", "")
        .replace("Mići princ", "Mići Princ"),
    )

for event in exb.findall(".//event"):
    try:
        event.text = event.text.strip()
    except:
        # print(event.text, type(event.text))
        event.getparent().remove(event)

timeline = {i.get("id"): float(i.get("time")) for i in exb.findall(".//tli")}

events = sorted(
    [i for i in exb.findall(".//event") if i.text.strip() != "*"],
    key=lambda event: timeline[event.get("start")],
)


text_in_exb = " ".join(
    [i.text for i in events]
)  # .replace("  ", " ").replace("  ", " ")
text_in_txt = Path(transcript).read_text()  # .replace("  ", " ").replace("  ", " ")
try:
    assert text_in_exb == text_in_txt, "Texts in txt and EXB differ!"
except AssertionError:
    for i in range(max(len(text_in_exb), len(text_in_txt))):
        if text_in_exb[:i] == text_in_txt[:i]:
            continue
        else:
            print(
                f"""
                  In txt: {text_in_txt[i-30:i+30]}
                  In exb: {text_in_exb[i-30:i+30]}
                  
                  """
            )
            break
    exit()
events_to_delete = list()
for event in exb.findall(".//event"):
    if event.text.strip() != "*":
        continue
    # print("Found an asterisk!", ET.tostring(event).decode("utf8"))
    candidate_events = exb.findall(f".//event[@start='{event.get('start')}'][@end='{event.get('end')}']")
    event_to_move = [e for e in candidate_events if e != event][0]
    event.text = event_to_move.text
    
    events_to_delete.append(event_to_move)
for event in events_to_delete:
    event.getparent().remove(event)
    
tiers_to_delete = []
for tier in exb.findall(".//tier"):
    if len(tier.findall(".//event")) == 0:
        tiers_to_delete.append(tier)
for tier in tiers_to_delete:
    tier.getparent().remove(tier)

print("Saving to ", outfile)
ET.indent(exb, space="\t")
exb.getroottree().write(
    Path(outfile),
    pretty_print=True,
    encoding="utf8",
    xml_declaration='<?xml version="1.0" encoding="UTF-8"?>',
)

Path(outfile).write_text(
    Path(outfile)
    .read_text()
    .replace(
        "<?xml version='1.0' encoding='UTF8'?>",
        '<?xml version="1.0" encoding="UTF-8"?>',
    )
)
