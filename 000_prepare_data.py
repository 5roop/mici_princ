from pathlib import Path
import string
from utils import fix_transcription, fix_audio, prep_helper_files, SAMPLEID

input_dir = Path("input_data")
kaldi_dir = Path("kaldi_data")
kaldi_outdir = Path("./kaldi_output")
kaldi_outdir.mkdir(exist_ok=True)
models = Path("./models")
assert models.exists(), "No folder named models in CWD"
kaldi = Path("/opt/kaldi")
assert kaldi.exists(), f"Kaldi not found at {kaldi}"
assert Path("./lexicon.py").exists(), f"Lexicon.py not found in CWD"
python = Path("/home/rupnik/anaconda3/envs/p37/bin/python")
for file in list(input_dir.glob("*.wav")):
    transcript_inpath = input_dir / file.with_suffix(".txt").name
    transcript_outpath = kaldi_dir / file.with_suffix(".txt").name
    audio_inpath = input_dir / file.with_suffix(".wav").name
    audio_outpath = kaldi_dir / file.with_suffix(".wav").name
    spk2utt = kaldi_dir / file.with_suffix(".spk2utt").name
    wavscp = kaldi_dir / file.with_suffix(".wavscp").name
    ali = kaldi_dir / file.with_suffix(".ctm").name
    fix_transcription(transcript_inpath, transcript_outpath)
    fix_audio(audio_inpath, audio_outpath)
    prep_helper_files(wavscp, spk2utt, audio_outpath)

    from subprocess import run

    cmd = f"""
    wavscp={wavscp}
    spk2utt={spk2utt}
    text={transcript_outpath}
    tmpdir={kaldi_outdir}

    models={models}
    kaldi={kaldi}

    export LC_ALL=C

    mkdir -p $tmpdir

    cut -f2- -d' ' $text | tr ' ' '\n' | sort -u > $tmpdir/wlist
    {python} lexicon.py $tmpdir/wlist $models/phonetisaurus-hr/model.fst $tmpdir
    $kaldi/src/featbin/compute-mfcc-feats  --config=$models/nnet3/conf/mfcc.conf scp:$wavscp ark:$tmpdir/mfcc.ark 
    $kaldi/src/online2bin/ivector-extract-online2 --config=$models/nnet3/conf/ivector.conf ark:$spk2utt ark:$tmpdir/mfcc.ark ark:$tmpdir/ivec.ark
    $kaldi/egs/wsj/s5/utils/sym2int.pl -f 2- $tmpdir/words.txt $text > $tmpdir/text.int
    $kaldi/src/bin/compile-train-graphs $models/nnet3/tdnn1a_sp/tree $models/nnet3/tdnn1a_sp/final.mdl $tmpdir/L.fst ark:$tmpdir/text.int ark:$tmpdir/graphs.fsts
    $kaldi/src/nnet3bin/nnet3-latgen-faster --allow_partial=True --online-ivectors=ark:$tmpdir/ivec.ark --online-ivector-period=10 $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/graphs.fsts ark:$tmpdir/mfcc.ark ark:$tmpdir/ali.lat
    $kaldi/src/latbin/lattice-align-words $tmpdir/word_boundary.int $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:- | $kaldi/src/latbin/lattice-to-ctm-conf ark:- - | $kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/words.txt - > $tmpdir/ali.ctm
    cp $tmpdir/ali.ctm {ali}

    # # # My additions:
    $kaldi/src/latbin/lattice-align-phones --replace-output-symbols=true $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:$tmpdir/phone_aligned.lats
    $kaldi/src/latbin/lattice-to-ctm-conf --inv-acoustic-scale=10 --decode-mbr ark:$tmpdir/phone_aligned.lats $tmpdir/trans_phones.ctm
    """

    run(cmd, shell=True)
