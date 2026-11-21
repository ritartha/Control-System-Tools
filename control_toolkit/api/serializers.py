"""DRF Serializers for all API endpoints."""
from rest_framework import serializers


class TransferFunctionSerializer(serializers.Serializer):
    num = serializers.CharField()
    den = serializers.CharField()
    unity = serializers.BooleanField(default=True)


class TimeResponseSerializer(serializers.Serializer):
    num = serializers.CharField()
    den = serializers.CharField()
    unity = serializers.BooleanField(default=True)
    signal_type = serializers.ChoiceField(choices=["Step", "Ramp", "Impulse"], default="Step")
    t_final = serializers.FloatField(required=False, allow_null=True, default=None)


class FrequencyResponseSerializer(serializers.Serializer):
    num = serializers.CharField()
    den = serializers.CharField()
    unity = serializers.BooleanField(default=True)


class RootLocusSerializer(serializers.Serializer):
    num = serializers.CharField()
    den = serializers.CharField()
    k_max = serializers.FloatField(required=False, allow_null=True, default=None)
    k_eval = serializers.FloatField(required=False, allow_null=True, default=None)
    zeta_target = serializers.FloatField(required=False, allow_null=True, default=None)


class RouthHurwitzSerializer(serializers.Serializer):
    num = serializers.CharField()
    den = serializers.CharField()
    unity = serializers.BooleanField(default=True)
    K = serializers.FloatField(default=1.0)


class NicholsSerializer(serializers.Serializer):
    num = serializers.CharField()
    den = serializers.CharField()
    target_mr = serializers.FloatField(required=False, allow_null=True, default=None)


class PIDSerializer(serializers.Serializer):
    num = serializers.CharField()
    den = serializers.CharField()
    Kp = serializers.FloatField(default=1.0)
    Ki = serializers.FloatField(default=0.0)
    Kd = serializers.FloatField(default=0.0)
    N = serializers.FloatField(default=20.0)
    zn_type = serializers.ChoiceField(
        choices=["P", "PI", "PID"], required=False, allow_null=True, default=None
    )


class CompensatorSerializer(serializers.Serializer):
    num = serializers.CharField()
    den = serializers.CharField()
    comp_type = serializers.ChoiceField(choices=["Lead", "Lag", "Lag-Lead"], default="Lead")
    K = serializers.FloatField(default=1.0)
    pm = serializers.FloatField(default=45.0)
    margin = serializers.FloatField(default=5.0)
    w1 = serializers.FloatField(required=False, allow_null=True, default=None)
