(function() {
    'use strict';

    function resultRow(label, value) {
        return '<div class="result-row"><span class="label">' + label + '</span><span class="value">' + value + '</span></div>';
    }

    function buildSpecsHTML(label, specs) {
        if (!specs) return '';
        var html = '<div class="result-section">';
        html += '<h4>' + label + '</h4>';
        html += resultRow('Stable', specs.stable ? '<span class="stable">YES ✓</span>' : '<span class="unstable">NO ✗</span>');
        if (specs.final_value !== undefined) html += resultRow('Final Value', App.formatNumber(specs.final_value));
        if (specs.ess !== undefined) html += resultRow('eₛₛ', App.formatNumber(specs.ess));
        if (specs.tr !== undefined) html += resultRow('Rise Time (tᵣ)', App.formatNumber(specs.tr) + ' s');
        if (specs.tp !== undefined) html += resultRow('Peak Time (tₚ)', App.formatNumber(specs.tp) + ' s');
        if (specs.Mp !== undefined) html += resultRow('Overshoot (Mₚ)', App.formatNumber(specs.Mp) + ' %');
        if (specs.ts !== undefined) html += resultRow('Settling Time (tₛ)', App.formatNumber(specs.ts) + ' s');
        if (specs.poles && Array.isArray(specs.poles)) {
            html += '<h4 style="margin-top:8px;">Closed-Loop Poles</h4>';
            specs.poles.forEach(function(p) { html += resultRow('', App.formatValue(p)); });
        }
        html += '</div>';
        return html;
    }

    function buildResultHTML(data) {
        var r = data.results;
        var html = '';

        // Z-N info
        if (r.zn_info) {
            html += '<div class="result-section">';
            html += '<h4>Ziegler-Nichols Info</h4>';
            if (typeof r.zn_info === 'string') {
                html += resultRow('', r.zn_info);
            } else if (typeof r.zn_info === 'object') {
                for (var key in r.zn_info) {
                    if (r.zn_info.hasOwnProperty(key)) {
                        html += resultRow(key, App.formatNumber(r.zn_info[key]));
                    }
                }
            }
            html += '</div>';
        }

        // PID gains
        html += '<div class="result-section">';
        html += '<h4>PID Controller</h4>';
        html += resultRow('Kp', App.formatNumber(r.Kp));
        html += resultRow('Ki', App.formatNumber(r.Ki));
        html += resultRow('Kd', App.formatNumber(r.Kd));
        html += '</div>';

        // Before/After comparison
        html += buildSpecsHTML('Without PID (Original)', r.without_pid);
        html += buildSpecsHTML('With PID (Controlled)', r.with_pid);

        return html;
    }

    function calculate(znType) {
        var system = App.getSystem();
        var kpInput = document.getElementById('pid-kp');
        var kiInput = document.getElementById('pid-ki');
        var kdInput = document.getElementById('pid-kd');
        var nInput = document.getElementById('pid-n');

        var params = {
            num: system.num,
            den: system.den,
            unity: system.unity
        };

        if (znType) {
            params.zn_type = znType;
        }

        var kp = kpInput ? parseFloat(kpInput.value) : NaN;
        var ki = kiInput ? parseFloat(kiInput.value) : NaN;
        var kd = kdInput ? parseFloat(kdInput.value) : NaN;
        var N = nInput ? parseFloat(nInput.value) : NaN;

        if (!isNaN(kp)) params.Kp = kp;
        if (!isNaN(ki)) params.Ki = ki;
        if (!isNaN(kd)) params.Kd = kd;
        if (!isNaN(N)) params.N = N;

        var infoPanel = document.getElementById('pid-info');
        App.showLoading('pid-info');

        App.apiCall('pid-tuning', params).then(function(data) {
            App.hideLoading('pid-info');
            if (data.plots && data.plots.step_response) {
                Plots.renderPlot('pid-plot', data.plots.step_response);
            }

            var r = data.results;
            // Update PID gain inputs if auto-tuned
            if (r) {
                if (r.Kp !== undefined && kpInput) kpInput.value = r.Kp;
                if (r.Ki !== undefined && kiInput) kiInput.value = r.Ki;
                if (r.Kd !== undefined && kdInput) kdInput.value = r.Kd;
            }

            if (r && infoPanel) {
                infoPanel.innerHTML = buildResultHTML(data);
            }
        }).catch(function(err) {
            App.hideLoading('pid-info');
            App.showToast(err.message);
        });
    }

    // Manual simulate button
    var calcBtn = document.getElementById('pid-calc-btn');
    if (calcBtn) {
        calcBtn.addEventListener('click', function() {
            calculate(null);
        });
    }

    // Ziegler-Nichols auto-tune buttons
    var znBtns = document.querySelectorAll('.zn-btn');
    znBtns.forEach(function(btn) {
        btn.addEventListener('click', function() {
            var znType = btn.getAttribute('data-zn');
            calculate(znType);
        });
    });
})();
