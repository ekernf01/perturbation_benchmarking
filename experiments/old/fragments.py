

cellFateClassifier = sklearn.linear_model.LogisticRegressionCV(Cs=10,
                                                               solver="liblinear",
                                                               verbose=1)
cellFateClassifier.fit(X=ko_lab_esc_data.X, y=ko_lab_esc_data.obs["leiden"])
import joblib
joblib.dump(cellFateClassifier, os.environ["PERTURBATION_PATH"] + "/nakatake/" + "classifier.joblib")
cellFateClassifier = joblib.load(os.environ["PERTURBATION_PATH"] + "/nakatake/" + "classifier.joblib") 

