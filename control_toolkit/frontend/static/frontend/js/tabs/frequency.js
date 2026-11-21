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
        html += resultRow('Gain Crossover (ωgc)', App.formatNumber(r.wgc) + ' rad/s');
        html += resultRow('Phase Crossover (ωpc)', App.formatNumber(r.wpc) + ' rad/s');
        html += '</div>';

        // Stability verdict
        html += '<div class="result-section">';
        html += '<h4>Stability</h4>';
        html += resultRow('Closed-Loop Stable',
            r.closed_loop_stable ? '<span class="stable">YES ✓</span>' : '<span class="unstable">NO ✗</span>');
        if (r.margins_verdict) {
            html += resultRow('Verdict', r.margins_verdict);
        }
        html += '</div>';

        // Closed-loop performance
        html += '<div class="result-section">';
        html += '<h4>Closed-Loop Performance</h4>';
        html += resultRow('Bandwidth (ωBW)', App.formatNumber(r.bandwidth) + ' rad/s');
        html += resultRow('DC Gain', App.formatNumber(r.dc_gain));
        if (r.Mr !== undefined && r.Mr !== null) {
            html += resultRow('Resonant Peak (Mr)', App.formatNumber(r.Mr) + ' dB');
        }
        if (r.wr !== undefined && r.wr !== null) {
            html += resultRow('Resonant Freq (ωr)', App.formatNumber(r.wr) + ' rad/s');
        }
        if (r.equiv_zeta !== undefined && r.equiv_zeta !== null) {
            html += resultRow('Equivalent ζ', App.formatNumber(r.equiv_zeta));
        }
        html += '</div>';

        return html;
    }

    function calculate() {
        var system = App.getSystem();
        var params = {
            num: system.num,
            den: system.den,
            unity: system.unity
        };

        var infoPanel = document.getElementById('freq-info');
        App.showLoading('freq-info');

        App.apiCall('frequency-response', params).then(function(data) {
            App.hideLoading('freq-info');
            if (data.plots) {
                if (data.plots.bode_mag) Plots.renderPlot('freq-bode-mag', data.plots.bode_mag);
                if (data.plots.bode_phase) Plots.renderPlot('freq-bode-phase', data.plots.bode_phase);
                if (data.plots.nyquist) Plots.renderPlot('freq-nyquist', data.plots.nyquist);
            }
            if (data.results && infoPanel) {
                infoPanel.innerHTML = buildResultHTML(data);
            }
        }).catch(function(err) {
            App.hideLoading('freq-info');
            App.showToast(err.message);
        });
    }

    var calcBtn = document.getElementById('freq-calc-btn');
    if (calcBtn) calcBtn.addEventListener('click', calculate);
})();
