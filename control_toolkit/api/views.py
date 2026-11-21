"""DRF API Views - each endpoint validates input and delegates to engine."""
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from . import engine
from .serializers import (
    TransferFunctionSerializer, TimeResponseSerializer,
    FrequencyResponseSerializer, RootLocusSerializer,
    RouthHurwitzSerializer, NicholsSerializer,
    PIDSerializer, CompensatorSerializer,
)


class TransferFunctionView(APIView):
    def post(self, request):
        ser = TransferFunctionSerializer(data=request.data)
        ser.is_valid(raise_exception=True)
        try:
            result = engine.compute_transfer_function(
                ser.validated_data["num"],
                ser.validated_data["den"],
                ser.validated_data["unity"],
            )
            return Response(result)
        except Exception as e:
            return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)


class TimeResponseView(APIView):
    def post(self, request):
        ser = TimeResponseSerializer(data=request.data)
        ser.is_valid(raise_exception=True)
        try:
            result = engine.compute_time_response(
                ser.validated_data["num"],
                ser.validated_data["den"],
                ser.validated_data["unity"],
                ser.validated_data["signal_type"],
                ser.validated_data["t_final"],
            )
            return Response(result)
        except Exception as e:
            return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)


class FrequencyResponseView(APIView):
    def post(self, request):
        ser = FrequencyResponseSerializer(data=request.data)
        ser.is_valid(raise_exception=True)
        try:
            result = engine.compute_frequency_response(
                ser.validated_data["num"],
                ser.validated_data["den"],
                ser.validated_data["unity"],
            )
            return Response(result)
        except Exception as e:
            return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)


class RootLocusView(APIView):
    def post(self, request):
        ser = RootLocusSerializer(data=request.data)
        ser.is_valid(raise_exception=True)
        try:
            result = engine.compute_root_locus(
                ser.validated_data["num"],
                ser.validated_data["den"],
                ser.validated_data["k_max"],
                ser.validated_data["k_eval"],
                ser.validated_data["zeta_target"],
            )
            return Response(result)
        except Exception as e:
            return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)


class RouthHurwitzView(APIView):
    def post(self, request):
        ser = RouthHurwitzSerializer(data=request.data)
        ser.is_valid(raise_exception=True)
        try:
            result = engine.compute_routh_hurwitz(
                ser.validated_data["num"],
                ser.validated_data["den"],
                ser.validated_data["unity"],
                ser.validated_data["K"],
            )
            return Response(result)
        except Exception as e:
            return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)


class NicholsView(APIView):
    def post(self, request):
        ser = NicholsSerializer(data=request.data)
        ser.is_valid(raise_exception=True)
        try:
            result = engine.compute_nichols(
                ser.validated_data["num"],
                ser.validated_data["den"],
                ser.validated_data["target_mr"],
            )
            return Response(result)
        except Exception as e:
            return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)


class PIDView(APIView):
    def post(self, request):
        ser = PIDSerializer(data=request.data)
        ser.is_valid(raise_exception=True)
        try:
            result = engine.compute_pid(
                ser.validated_data["num"],
                ser.validated_data["den"],
                ser.validated_data["Kp"],
                ser.validated_data["Ki"],
                ser.validated_data["Kd"],
                ser.validated_data["N"],
                ser.validated_data["zn_type"],
            )
            return Response(result)
        except Exception as e:
            return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)


class CompensatorView(APIView):
    def post(self, request):
        ser = CompensatorSerializer(data=request.data)
        ser.is_valid(raise_exception=True)
        try:
            result = engine.compute_compensator(
                ser.validated_data["num"],
                ser.validated_data["den"],
                ser.validated_data["comp_type"],
                ser.validated_data["K"],
                ser.validated_data["pm"],
                ser.validated_data["margin"],
                ser.validated_data["w1"],
            )
            return Response(result)
        except Exception as e:
            return Response({"error": str(e)}, status=status.HTTP_400_BAD_REQUEST)
