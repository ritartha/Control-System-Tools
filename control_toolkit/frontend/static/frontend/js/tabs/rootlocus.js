(function() {
    'use strict';

    function resultRow(label, value) {
        return '<div class="result-row"><span class="label">' + label + '</span><span class="value">' + value + '</span></div>';
    }

    function buildResultHTML(data) {
        var r = data.results;
        var html = '';

        // Open-loop poles and zeros
        html += '<div class="result-section">';
        html += '<h4>Open-Loop Poles</h4>';
        if (r.poles && r.poles.length) {
            r.poles.forEach(function(p) { html += resultRow('', App.formatValue(p)); });
        } else {
            html += resultRow('', 'None');
        }
        html += '</div>';

        html += '<div class="result-section">';
        html += '<h4>Open-Loop Zeros</h4>';
        if (r.zeros && r.zeros.length) {
            r.zeros.forEach(function(z) { html += resultRow('', App.formatValue(z)); });
        } else {
            html += resultRow('', 'None');
        }
        html += '</div>';

        // Asymptotes
        html += '<div class="result-section">';
        html += '<h4>Asymptotes</h4>';
        html += resultRow('Branches', App.formatNumber(r.branches));
        html += resultRow('# Asymptotes', App.formatNumber(r.num_asymptotes));
        if (r.asymptote_angles && r.asymptote_angles.length) {
            html += resultRow('Angles', r.asymptote_angles.map(function(a) { return App.formatNumber(a) + '°'; }).join(', '));
        }
        if (r.centroid !== undefined && r.centroid !== null) {
            html += resultRow('Centroid (σₐ)', App.formatNumber(r.centroid));
        }
        html += '</div>';

        // Breakpoints
        if (r.breakpoints && r.breakpoints.length) {
            html += '<div class="result-section">';
            html += '<h4>Breakaway/Break-in Points</h4>';
            r.breakpoints.forEach(function(bp) {
                html += resultRow('', App.formatValue(bp));
            });
            html += '</div>';
        }

        // jω crossing
        if (r.jw_crossing) {
            html += '<div class="result-section">';
            html += '<h4>jω-Axis Crossing</h4>';
            if (typeof r.jw_crossing === 'object' && !Array.isArray(r.jw_crossing)) {
                if (r.jw_crossing.K !== undefined) html += resultRow('K', App.formatNumber(r.jw_crossing.K));
                if (r.jw_crossing.w !== undefined) html += resultRow('ω', App.formatNumber(r.jw_crossing.w) + ' rad/s');
            } else if (Array.isArray(r.jw_crossing)) {
                r.jw_crossing.forEach(function(c) { html += resultRow('', c); });
            } else {
                html += resultRow('', r.jw_crossing);
            }
            html += '</div>';
        }

        // Stable range
        if (r.stable_range) {
            html += '<div class="result-section">';
            html += '<h4>Stable K Range</h4>';
            html += resultRow('Range', r.stable_range);
            html += '</div>';
        }

        // Zeta target points
        if (r.zeta_points && r.zeta_points.length) {
            html += '<div class="result-section">';
            html += '<h4>ζ Target Points</h4>';
            r.zeta_points.forEach(function(zp) {
                html += resultRow('', App.formatValue(zp));
            });
            html += '</div>';
        }

        // K evaluation
        if (r.k_eval_poles) {
            html += '<div class="result-section">';
            html += '<h4>Poles at Evaluated K</h4>';
            if (Array.isArray(r.k_eval_poles)) {
                r.k_eval_poles.forEach(function(p) { html += resultRow('', App.formatValue(p)); });
            } else {
                html += resultRow('', r.k_eval_poles);
            }
            html += '</div>';
        }

        return html;
    }

    function calculate() {
        var system = App.getSystem();

        var kMaxInput = document.getElementById('rlocus-kmax');
        var kEvalInput = document.getElementById('rlocus-keval');
        var zetaInput = document.getElementById('rlocus-zeta');

        var params = {
            num: system.num,
            den: system.den,
            unity: system.unity
        };

        var kMax = kMaxInput ? parseFloat(kMaxInput.value) : NaN;
        if (!isNaN(kMax)) params.k_max = kMax;

        var kEval = kEvalInput ? parseFloat(kEvalInput.value) : NaN;
        if (!isNaN(kEval)) params.k_eval = kEval;

        var zetaTarget = zetaInput ? parseFloat(zetaInput.value) : NaN;
        if (!isNaN(zetaTarget)) params.zeta_target = zetaTarget;

        var infoPanel = document.getElementById('rlocus-info');
        App.showLoading('rlocus-info');

        App.apiCall('root-locus', params).then(function(data) {
            App.hideLoading('rlocus-info');
            if (data.plots && data.plots.root_locus) {
                Plots.renderPlot('rlocus-plot', data.plots.root_locus);
            }
            if (data.results && infoPanel) {
                infoPanel.innerHTML = buildResultHTML(data);
            }
        }).catch(function(err) {
            App.hideLoading('rlocus-info');
            App.showToast(err.message);
        });
    }

    var calcBtn = document.getElementById('rlocus-calc-btn');
    if (calcBtn) calcBtn.addEventListener('click', calculate);
})();
