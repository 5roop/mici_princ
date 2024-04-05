from pathlib import Path

available_fixed_exbs = Path("/home/nikola/transfer/mici_princ/").glob("MP*.exb")
expected_files = [Path("inspected_exbs_fixed")/ i.name for i in available_fixed_exbs]

rule GatherChecked:
    input: expected_files

rule Check:
    input:
        exb = "/home/nikola/transfer/mici_princ/{file}.exb",
        transcript = "/home/nikola/transfer/mici_princ/{file}.txt.txt",
        jsonfile = "/home/nikola/transfer/mici_princ/{file}.txt.json",
    output:
        out = "inspected_exbs_fixed/{file}.exb"
    script:
        "sanitycheck.py"

rule GatherJsons:
    default_target: True
    input:
        [Path("aligned_jsons")/ i.with_suffix(".json").name for i in expected_files]
    run:
        import json
        from pathlib import Path
        from colorama import Fore, Back, Style
        speakers = set()
        speakers_in_files = dict()
        for file in input:
            f = json.loads(Path(file).read_text())
            for u in f:
                speaker = u["speaker"]
                speakers = {speaker, *speakers}
                fname = Path(file).name
                speakers_in_files[speaker] = {fname, *speakers_in_files.get(speaker, set())}
        print(Fore.YELLOW, speakers, end="\n\n\n")
        for s, f in speakers_in_files.items():
            print(Fore.YELLOW, s, Fore.CYAN, f)
    
rule DoJSON:
    input:
        exb="inspected_exbs_fixed/{file}.exb",
        transcript="/home/nikola/transfer/mici_princ/{file}.txt.txt",
        gold="/home/nikola/transfer/mici_princ/{file}.txt",
        jsonfile="/home/nikola/transfer/mici_princ/{file}.txt.json",
    output:
        out="aligned_jsons/{file}.json"
    script:
        "generate_jsons.py"