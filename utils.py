from pathlib import Path
import string

SAMPLEID = "sampleid"


def fix_transcription(transcript_inpath: Path, transcript_outpath: Path) -> None:
    t = transcript_inpath.read_text()

    # Delete newlines:
    t = " ".join(t.split())

    t = t.replace("ȋ", "i").replace("ȅ", "e")
    # Fix numbers:
    words = t.split()
    for i, word in enumerate(words):
        contains_digits = any([i in string.digits for i in word])
        if not contains_digits:
            continue
        context = 4
        start_index = max([0, i - context])
        end_index = min([len(words), i + context])
        print(" ".join(words[start_index:end_index]))
        correction = input("Correction?: ")
        words[i] = correction
    t = " ".join(words)

    # Fix punctuation:
    for l in string.punctuation + "„“•–":
        t = t.replace(l, "")
    # Fix capital letters:
    t = t.casefold()
    t = " ".join(t.split())
    print("In finished text we have the following characters:", sorted(list(set(t))))

    transcript_outpath.write_text(SAMPLEID + " " + t)


def fix_audio(audio_inpath: Path, audio_outpath: Path) -> None:
    import subprocess

    trim = "" if audio_inpath.name != "MP_02.wav" else "-ss 62.4"
    cmd = f"""ffmpeg {trim} -i {audio_inpath} -ac 1 -ar 16000 {audio_outpath}"""
    print("\n\n\n\n", cmd)
    subprocess.run(
        cmd,
        shell=True,
    )


def prep_helper_files(wavscp: Path, spk2utt: Path, audio_outpath) -> None:
    wavscp.write_text(SAMPLEID + " " + str(audio_outpath))
    spk2utt.write_text("speaker " + SAMPLEID)
