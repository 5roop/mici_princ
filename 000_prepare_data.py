from pathlib import Path
from tqdm import tqdm
import string
from utils import fix_transcription, fix_audio, prep_helper_files, SAMPLEID
from subprocess import run

input_data = Path("/home/nikola/transfer/mici_princ")
downsampled = Path("downsampled_and_trimmed_wavs")
MPs = sorted(list(input_data.glob("*.txt.txt")))

for MP in tqdm(MPs):
    MP = MP.name.replace(".txt.txt", "")
    print(MP)
    transcript = input_data / f"{MP}.txt.txt"
    wav = downsampled / f"{MP}.wav"

    (Path("kaldi_data") / "text").write_text(
        f"""{SAMPLEID} """ + transcript.read_text()
    )
    (Path("kaldi_data") / "wav.scp").write_text(
        f"""{SAMPLEID} {downsampled /( MP+'.wav')}"""
    )

    run("bash 001_kaldi.sh", shell=True, capture_output=True)
    ali = Path("kaldi_output/ali.ctm").read_text()
    if len(ali.rstrip().split("\n")) > 1:
        # print("Worked for", MP)
        Path(f"kaldi_output/{MP}.ali.ctm").write_text(ali)
    else:
        print("Failed for", MP)
        continue
