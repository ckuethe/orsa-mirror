include(../../../orsa.pri)

TEMPLATE = subdirs

SUBDIRS += SurveyReview.app.pro 

unix:!macx {
	SUBDIRS += ExtractObservations.app.pro
	SUBDIRS += ProcessObservations.app.pro
	SUBDIRS += GenMake.app.pro
	SUBDIRS += GenSkipFile.app.pro
	SUBDIRS += DetectionEfficiency.app.pro
	SUBDIRS += DetectionEfficiencyFit.app.pro
	SUBDIRS += DetectionEfficiencyPlot.app.pro
	SUBDIRS += SurveyReviewWorkGenerator.app.pro
	SUBDIRS += SurveyReviewMultipleJobsSubmission.app.pro
#	SUBDIRS += SurveyReviewValidator.app.pro
#	SUBDIRS += SurveyReviewAssimilator.app.pro
	SUBDIRS += SurveyReviewMerge.app.pro
	SUBDIRS += SurveyReviewMergeNEOs.app.pro
	SUBDIRS += InspectMerge.app.pro
	SUBDIRS += InspectMergePlot.app.pro
	SUBDIRS += InspectMergePlot_ae.app.pro
	SUBDIRS += InspectMergePlot_ai.app.pro
	SUBDIRS += InspectMergePlot_aL.app.pro
	SUBDIRS += InspectMergePlot_sky.app.pro
	SUBDIRS += InspectMergePlot_xy.app.pro
	SUBDIRS += Inspect.app.pro
	SUBDIRS += InspectPlot.app.pro
	SUBDIRS += InspectPlot2D.app.pro
	SUBDIRS += InspectPlot2DLunar.app.pro
	SUBDIRS += LunationCopy.app.pro
}
