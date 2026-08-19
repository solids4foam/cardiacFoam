from pathlib import Path
from openfoam_driver.plugins.cardiacfoam.artifacts_predictor import predict_cardiac_artifacts

arts = predict_cardiac_artifacts(Path("../../tutorials/NiedererEtAl2011/NiedererEtAl2011verification"), None)
print("Predicted artifacts for ../../tutorials:")
for a in arts:
    print(f"- {a.artifact_id}")

arts = predict_cardiac_artifacts(Path("../../../tutorials/NiedererEtAl2011/NiedererEtAl2011verification"), None)
print("Predicted artifacts for ../../../tutorials:")
for a in arts:
    print(f"- {a.artifact_id}")
