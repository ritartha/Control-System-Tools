(function() {
    'use strict';

    function resultRow(label, value) {
        return '<div class="result-row"><span class="label">' + label + '</span><span class="value">' + value + '</span></div>';
    }

    function buildRouthTable(routhTable, rowLabels) {
        if (!routhTable || !routhTable.length) return '';

        var html = '<table class="routh-table">';
        for (var i = 0; i < routhTable.length; i++) {
            html += '<tr>';
            if (rowLabels && rowLabels[i]) {
                html += '<th>' + rowLabels[i] + '</th>';
            }
            var row = routhTable[i];
            for (var j = 0; j < row.length; j++) {
                html += '<td>' + App.formatNumber(row[j]) + '</td>';
            }
            html += '</tr>';
        }
        html += '</table>';
        return html;
    }

    function buildResultHTML(data) {
        var r = data.results;
        var html = '';

        // Equation description
        if (r.equation_desc) {
            html += '<div class="result-section">';
            html += '<h4>Characteristic Equation</h4>';
            html += resultRow('', r.equation_desc);
            html += '</div>';
        }

        // Routh table
        html += '<div class="result-section">';
        html += '<h4>Routh Array</h4>';
        html += buildRouthTable(r.routh_table, r.row_labels);
        html += '</div>';

        // Verdict
        html += '<div class="result-section">';
        html += '<h4>Verdict</h4>';
        if (r.verdict) {
            var verdictClass = r.verdict.toLowerCase().indexOf('stable') !== -1 &&
                               r.verdict.toLowerCase().indexOf('unstable') === -1 ? 'stable' : 'unstable';
            html += resultRow('Result', '<span class="' + verdictClass + '">' + r.verdict + '</span>');
        }
        if (r.sign_changes !== undefined) {
            html += resultRow('Sign Changes', App.formatNumber(r.sign_changes));
        }
        if (r.rhp_roots !== undefined) {
            html += resultRow('RHP Roots', App.formatNumber(r.rhp_roots));
        }
        if (r.jw_roots && r.jw_roots.length) {
            html += resultRow('jω-Axis Roots', r.jw_roots.join(', '));
        }
        html += '</div>';

        // Stable ranges
        if (r.stable_ranges) {
            html += '<div class="result-section">';
            html += '<h4>Stable Ranges</h4>';
            if (Array.isArray(r.stable_ranges)) {
                r.stable_ranges.forEach(function(sr) {
                    html += resultRow('', App.formatValue(sr));
                });
            } else {
                html += resultRow('', r.stable_ranges);
            }
            html += '</div>';
        }

        // Actual roots
        if (r.actual_roots && r.actual_roots.length) {
            html += '<div class="result-section">';
            html += '<h4>Actual Roots</h4>';
            r.actual_roots.forEach(function(root) {
                html += resultRow('', App.formatValue(root));
            });
            html += '</div>';
        }

        // Notes
        if (r.notes && r.notes.length) {
            html += '<div class="result-section">';
            html += '<h4>Notes</h4>';
            r.notes.forEach(function(note) {
                html += resultRow('', note);
            });
            html += '</div>';
        }

        return html;
    }

    function calculate() {
        var system = App.getSystem();
        var kInput = document.getElementById('routh-K');

        var params = {
            num: system.num,
            den: system.den,
            unity: system.unity
        };

        var K = kInput ? kInput.value.trim() : '';
        if (K) params.K = parseFloat(K) || K;

        var infoPanel = document.getElementById('routh-info');
        App.showLoading('routh-info');

        App.apiCall('routh-hurwitz', params).then(function(data) {
            App.hideLoading('routh-info');
            if (data.plots && data.plots.roots) {
                Plots.renderPlot('routh-plot', data.plots.roots);
            }
            if (data.results && infoPanel) {
                infoPanel.innerHTML = buildResultHTML(data);
            }
        }).catch(function(err) {
            App.hideLoading('routh-info');
            App.showToast(err.message);
        });
    }

    var calcBtn = document.getElementById('routh-calc-btn');
    if (calcBtn) calcBtn.addEventListener('click', calculate);
})();
