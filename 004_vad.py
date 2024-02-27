from pyannote.audio import Model

model = Model.from_pretrained(
    "pyannote/segmentation", use_auth_token="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"
)
from pyannote.audio.pipelines import VoiceActivityDetection

pipeline = VoiceActivityDetection(segmentation=model)
HYPER_PARAMETERS = {
    # onset/offset activation thresholds
    "onset": 0.9,
    "offset": 0.5,
    # remove speech regions shorter than that many seconds.
    "min_duration_on": 0.0,
    # fill non-speech regions shorter than that many seconds.
    "min_duration_off": 0.0,
}
pipeline.instantiate(HYPER_PARAMETERS)
output = pipeline("input_data/MP_02.wav")

print(poo)
