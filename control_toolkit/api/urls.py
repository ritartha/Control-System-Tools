"""API URL routing."""
from django.urls import path
from . import views

urlpatterns = [
    path("transfer-function/", views.TransferFunctionView.as_view(), name="transfer-function"),
    path("time-response/", views.TimeResponseView.as_view(), name="time-response"),
    path("frequency-response/", views.FrequencyResponseView.as_view(), name="frequency-response"),
    path("root-locus/", views.RootLocusView.as_view(), name="root-locus"),
    path("routh-hurwitz/", views.RouthHurwitzView.as_view(), name="routh-hurwitz"),
    path("nichols/", views.NicholsView.as_view(), name="nichols"),
    path("pid-tuning/", views.PIDView.as_view(), name="pid-tuning"),
    path("compensator/", views.CompensatorView.as_view(), name="compensator"),
]
