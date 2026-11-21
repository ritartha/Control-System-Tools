(function() {
    'use strict';

    function resultRow(label, value) {
        return '<div class="result-row"><span class="label">' + label + '</span><span class="value">' + value + '</span></div>';
    }

    function buildSpecsHTML(label, specs) {
        if (!specs) return '';
        var html = '';
        html += '<h4>' + label + '</h4>';
        html += resultRow('PM', App.formatNumber(specs.PM) + '°');
        html += resultRow('GM', App.formatNumber(specs.GM) + ' dB');
        html += resultRow('ωgc', App.formatNumber(specs.wgc) + ' rad/s');
        html += resultRow('Stable', specs.stable ? '<span class="stable">YES ✓</span>' : '<span class="unstable">NO ✗</span>');
        html += resultRow('Final Value', App.formatNumber(specs.final_value));
        html += resultRow('eₛₛ', App.formatNumber(specs.ess));
        html += resultRow('Rise Time (tᵣ)', App.formatNumber(specs.tr) + ' s');
        html += resultRow('Peak Time (tₚ)', App.formatNumber(specs.tp) + ' s');
        html += resultRow('Overshoot (Mₚ)', App.formatNumber(specs.Mp) + ' %');
        html += resultRow('Settling Time (tₛ)', App.formatNumber(specs.ts) + ' s');
        return html;
    }

    function buildResultHTML(data) {
        var r = data.results;
        var html = '';

        // Compensator type and TF
        html += '<div class="result-section">';
        html += '<h4>Compensator</h4>';
        if (r.comp_type) html += resultRow('Type', r.comp_type);
        if (r.compensator_tf) {
            html += '<pre class="tf-fraction">' + r.compensator_tf + '</pre>';
        }
        html += '</div>';

        // Design info
        if (r.design_info) {
            html += '<div class="result-section">';
            html += '<h4>Design Information</h4>';
            if (typeof r.design_info === 'object' && !Array.isArray(r.design_info)) {
                for (var key in r.design_info) {
                    if (r.design_info.hasOwnProperty(key)) {
                        html += resultRow(key, App.formatNumber(r.design_info[key]));
                    }
                }
            } else if (Array.isArray(r.design_info)) {
                r.design_info.forEach(function(item) { html += resultRow('', item); });
            } else {
                html += resultRow('', r.design_info);
            }
            html += '</div>';
        }

        // Gain info
        if (r.gain_info) {
            html += '<div class="result-section">';
            html += '<h4>Gain Information</h4>';
            if (typeof r.gain_info === 'object' && !Array.isArray(r.gain_info)) {
                for (var key in r.gain_info) {
                    if (r.gain_info.hasOwnProperty(key)) {
                        html += resultRow(key, App.formatNumber(r.gain_info[key]));
                    }
                }
            } else {
                html += resultRow('', r.gain_info);
            }
            html += '</div>';
        }

        // Before / After comparison
        html += '<div class="result-section">';
        html += buildSpecsHTML('Before Compensation', r.before);
        html += buildSpecsHTML('After Compensation', r.after);
        html += '</div>';

        // Note
        if (r.note) {
            html += '<div class="result-section">';
            html += '<h4>Note</h4>';
            html += resultRow('', r.note);
            html += '</div>';
        }

        return html;
    }

    function calculate() {
        var system = App.getSystem();
        var compTypeSelect = document.getElementById('comp-type');
        var kInput = document.getElementById('comp-K');
        var pmInput = document.getElementById('comp-pm');
        var marginInput = document.getElementById('comp-margin');
        var w1Input = document.getElementById('comp-w1');

        var params = {
            num: system.num,
            den: system.den,
            unity: system.unity
        };

        if (compTypeSelect) params.comp_type = compTypeSelect.value;

        var K = kInput ? parseFloat(kInput.value) : NaN;
        if (!isNaN(K)) params.K = K;

        var pm = pmInput ? parseFloat(pmInput.value) : NaN;
        if (!isNaN(pm)) params.pm = pm;

        var margin = marginInput ? parseFloat(marginInput.value) : NaN;
        if (!isNaN(margin)) params.margin = margin;

        var w1 = w1Input ? parseFloat(w1Input.value) : NaN;
        if (!isNaN(w1)) params.w1 = w1;
        else params.w1 = null;

        var infoPanel = document.getElementById('comp-info');
        App.showLoading('comp-info');

        App.apiCall('compensator', params).then(function(data) {
            App.hideLoading('comp-info');
            if (data.plots) {
                if (data.plots.bode_mag) Plots.renderPlot('comp-bode-mag', data.plots.bode_mag);
                if (data.plots.bode_phase) Plots.renderPlot('comp-bode-phase', data.plots.bode_phase);
                if (data.plots.step_response) Plots.renderPlot('comp-step', data.plots.step_response);
            }
            if (data.results && infoPanel) {
                infoPanel.innerHTML = buildResultHTML(data);
            }
        }).catch(function(err) {
            App.hideLoading('comp-info');
            App.showToast(err.message);
        });
    }

    var calcBtn = document.getElementById('comp-calc-btn');
    if (calcBtn) calcBtn.addEventListener('click', calculate);
})();
