(function() {
    'use strict';

    function resultRow(label, value) {
        return '<div class="result-row"><span class="label">' + label + '</span><span class="value">' + value + '</span></div>';
    }

    function buildResultHTML(data) {
        var r = data.results;
        var html = '';

        // Poles
        html += '<div class="result-section">';
        html += '<h4>Poles</h4>';
        if (r.poles && r.poles.length) {
            r.poles.forEach(function(p) {
                var pstr = typeof p === 'object' ? App.formatNumber(p.real) + (p.imag >= 0 ? ' + ' : ' - ') + 'j' + App.formatNumber(Math.abs(p.imag)) : p;
                html += resultRow('', pstr);
            });
        } else {
            html += resultRow('', 'None');
        }
        html += '</div>';

        // Zeros
        html += '<div class="result-section">';
        html += '<h4>Zeros</h4>';
        if (r.zeros && r.zeros.length) {
            r.zeros.forEach(function(z) {
                var zstr = typeof z === 'object' ? App.formatNumber(z.real) + (z.imag >= 0 ? ' + ' : ' - ') + 'j' + App.formatNumber(Math.abs(z.imag)) : z;
                html += resultRow('', zstr);
            });
        } else {
            html += resultRow('', 'None');
        }
        html += '</div>';

        // Stability
        html += '<div class="result-section">';
        html += '<h4>Stability</h4>';
        html += resultRow('Stable', r.stable ? '<span class="stable">YES ✓</span>' : '<span class="unstable">NO ✗</span>');
        html += '</div>';

        // Step specifications
        if (r.step_specs) {
            html += '<div class="result-section">';
            html += '<h4>Step Response Specifications</h4>';
            if (r.step_specs.rise_time !== undefined && r.step_specs.rise_time !== null) {
                html += resultRow('Rise Time (tᵣ)', App.formatNumber(r.step_specs.rise_time) + ' s');
            }
            if (r.step_specs.peak_time !== undefined && r.step_specs.peak_time !== null) {
                html += resultRow('Peak Time (tₚ)', App.formatNumber(r.step_specs.peak_time) + ' s');
            }
            if (r.step_specs.overshoot !== undefined && r.step_specs.overshoot !== null) {
                html += resultRow('Overshoot (Mₚ)', App.formatNumber(r.step_specs.overshoot) + ' %');
            }
            if (r.step_specs.settling_time !== undefined && r.step_specs.settling_time !== null) {
                html += resultRow('Settling Time (tₛ)', App.formatNumber(r.step_specs.settling_time) + ' s');
            }
            html += '</div>';
        }

        // Second order formulas
        if (r.second_order) {
            html += '<div class="result-section">';
            html += '<h4>2nd Order Parameters</h4>';
            if (r.second_order.wn !== undefined) {
                html += resultRow('ωₙ', App.formatNumber(r.second_order.wn) + ' rad/s');
            }
            if (r.second_order.zeta !== undefined) {
                html += resultRow('ζ', App.formatNumber(r.second_order.zeta));
            }
            if (r.second_order.wd !== undefined) {
                html += resultRow('ωd', App.formatNumber(r.second_order.wd) + ' rad/s');
            }
            if (r.second_order.sigma !== undefined) {
                html += resultRow('σ', App.formatNumber(r.second_order.sigma));
            }
            html += '</div>';
        }

        // Steady state error
        if (r.steady_state) {
            html += '<div class="result-section">';
            html += '<h4>Steady-State Error</h4>';
            if (r.steady_state.type !== undefined) {
                html += resultRow('System Type', 'Type ' + r.steady_state.type);
            }
            if (r.steady_state.Kp !== undefined) {
                html += resultRow('Kp', App.formatNumber(r.steady_state.Kp));
            }
            if (r.steady_state.Kv !== undefined) {
                html += resultRow('Kv', App.formatNumber(r.steady_state.Kv));
            }
            if (r.steady_state.Ka !== undefined) {
                html += resultRow('Ka', App.formatNumber(r.steady_state.Ka));
            }
            if (r.steady_state.ess_step !== undefined) {
                html += resultRow('eₛₛ (step)', App.formatNumber(r.steady_state.ess_step));
            }
            if (r.steady_state.ess_ramp !== undefined) {
                html += resultRow('eₛₛ (ramp)', App.formatNumber(r.steady_state.ess_ramp));
            }
            if (r.steady_state.ess_parabolic !== undefined) {
                html += resultRow('eₛₛ (parabolic)', App.formatNumber(r.steady_state.ess_parabolic));
            }
            html += '</div>';
        }

        return html;
    }

    function calculate() {
        var system = App.getSystem();
        var signalSelect = document.getElementById('time-signal');
        var tFinalInput = document.getElementById('time-tfinal');

        var params = {
            num: system.num,
            den: system.den,
            unity: system.unity,
            signal_type: signalSelect ? signalSelect.value : 'step'
        };

        var tFinal = tFinalInput ? parseFloat(tFinalInput.value) : NaN;
        if (!isNaN(tFinal) && tFinal > 0) {
            params.t_final = tFinal;
        }

        var infoPanel = document.getElementById('time-info');
        App.showLoading('time-info');

        App.apiCall('time-response', params).then(function(data) {
            App.hideLoading('time-info');
            if (data.plots) {
                if (data.plots.response) Plots.renderPlot('time-plot', data.plots.response);
                if (data.plots.pole_zero) Plots.renderPlot('time-pzmap', data.plots.pole_zero);
            }
            if (data.results && infoPanel) {
                infoPanel.innerHTML = buildResultHTML(data);
            }
        }).catch(function(err) {
            App.hideLoading('time-info');
            App.showToast(err.message);
        });
    }

    function loadSecondOrder() {
        var zetaInput = document.getElementById('time-zeta');
        var wnInput = document.getElementById('time-wn');

        var zeta = zetaInput ? parseFloat(zetaInput.value) : NaN;
        var wn = wnInput ? parseFloat(wnInput.value) : NaN;

        if (isNaN(zeta) || isNaN(wn)) {
            App.showToast('Please enter valid ζ and ωₙ values');
            return;
        }

        var numInput = document.getElementById('num-input');
        var denInput = document.getElementById('den-input');
        var unityCheckbox = document.getElementById('unity-checkbox');

        var wn2 = wn * wn;
        if (numInput) numInput.value = wn2.toString();
        if (denInput) denInput.value = '1 ' + (2 * zeta * wn).toString() + ' ' + wn2.toString();
        if (unityCheckbox) unityCheckbox.checked = false;

        // Trigger TF update
        if (numInput) numInput.dispatchEvent(new Event('input'));

        calculate();
    }

    // Initialize
    var calcBtn = document.getElementById('time-calc-btn');
    var load2Btn = document.getElementById('time-load2-btn');

    if (calcBtn) calcBtn.addEventListener('click', calculate);
    if (load2Btn) load2Btn.addEventListener('click', loadSecondOrder);
})();
