(function() {
    'use strict';

    function resultRow(label, value) {
        return '<div class="result-row"><span class="label">' + label + '</span><span class="value">' + value + '</span></div>';
    }

    function buildResultHTML(data) {
        var r = data.results;
        var html = '';

        // Stability Margins
        html += '<div class="result-section">';
        html += '<h4>Stability Margins</h4>';
        html += resultRow('Gain Margin (GM)', App.formatNumber(r.GM) + ' dB');
        html += resultRow('Phase Margin (PM)', App.formatNumber(r.PM) + '°');
        html += '</div>';

        // Stability
        html += '<div class="result-section">';
        html += '<h4>Stability</h4>';
        html += resultRow('Stable', r.stable ? '<span class="stable">YES ✓</span>' : '<span class="unstable">NO ✗</span>');
        html += '</div>';

        // Peak info
        if (r.peak_info) {
            html += '<div class="result-section">';
            html += '<h4>Peak Closed-Loop Response</h4>';
            html += resultRow('', r.peak_info);
            html += '</div>';
        }

        // Bandwidth
        if (r.bandwidth !== undefined && r.bandwidth !== null) {
            html += '<div class="result-section">';
            html += '<h4>Bandwidth</h4>';
            html += resultRow('ωBW', App.formatNumber(r.bandwidth) + ' rad/s');
            html += '</div>';
        }

        // Gain adjustments
        if (r.gain_adjustments) {
            html += '<div class="result-section">';
            html += '<h4>Gain Adjustments</h4>';
            if (typeof r.gain_adjustments === 'object' && !Array.isArray(r.gain_adjustments)) {
                for (var key in r.gain_adjustments) {
                    if (r.gain_adjustments.hasOwnProperty(key)) {
                        html += resultRow(key, App.formatNumber(r.gain_adjustments[key]));
                    }
                }
            } else if (Array.isArray(r.gain_adjustments)) {
                r.gain_adjustments.forEach(function(ga) {
                    html += resultRow('', ga);
                });
            } else {
                html += resultRow('', r.gain_adjustments);
            }
            html += '</div>';
        }

        return html;
    }

    function calculate() {
        var system = App.getSystem();
        var mrInput = document.getElementById('nichols-mr');

        var params = {
            num: system.num,
            den: system.den,
            unity: system.unity
        };

        var targetMr = mrInput ? parseFloat(mrInput.value) : NaN;
        if (!isNaN(targetMr)) params.target_mr = targetMr;

        var infoPanel = document.getElementById('nichols-info');
        App.showLoading('nichols-info');

        App.apiCall('nichols', params).then(function(data) {
            App.hideLoading('nichols-info');
            if (data.plots && data.plots.nichols) {
                Plots.renderPlot('nichols-plot', data.plots.nichols);
            }
            if (data.results && infoPanel) {
                infoPanel.innerHTML = buildResultHTML(data);
            }
        }).catch(function(err) {
            App.hideLoading('nichols-info');
            App.showToast(err.message);
        });
    }

    var calcBtn = document.getElementById('nichols-calc-btn');
    if (calcBtn) calcBtn.addEventListener('click', calculate);
})();
