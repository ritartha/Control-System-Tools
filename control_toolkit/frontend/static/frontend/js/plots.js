(function() {
    'use strict';

    function getDefaultLayout() {
        return {
            font: {
                family: 'Inter, system-ui, sans-serif',
                color: '#1e293b',
                size: 12
            },
            paper_bgcolor: '#ffffff',
            plot_bgcolor: '#f8fafc',
            margin: { l: 60, r: 30, t: 40, b: 50 },
            autosize: true,
            xaxis: {
                gridcolor: 'rgba(148,163,184,0.2)',
                zerolinecolor: 'rgba(148,163,184,0.4)'
            },
            yaxis: {
                gridcolor: 'rgba(148,163,184,0.2)',
                zerolinecolor: 'rgba(148,163,184,0.4)'
            }
        };
    }

    function deepMerge(target, source) {
        var result = {};
        var key;
        for (key in target) {
            if (target.hasOwnProperty(key)) {
                result[key] = target[key];
            }
        }
        for (key in source) {
            if (source.hasOwnProperty(key)) {
                if (
                    result[key] && typeof result[key] === 'object' && !Array.isArray(result[key]) &&
                    source[key] && typeof source[key] === 'object' && !Array.isArray(source[key])
                ) {
                    result[key] = deepMerge(result[key], source[key]);
                } else {
                    result[key] = source[key];
                }
            }
        }
        return result;
    }

    function renderPlot(elementId, plotData) {
        var el = document.getElementById(elementId);
        if (!el) return;

        var defaults = getDefaultLayout();
        var mergedLayout = deepMerge(defaults, plotData.layout || {});

        // Compute size from container, accounting for padding
        var rect = el.getBoundingClientRect();
        var style = window.getComputedStyle(el);
        var padX = parseFloat(style.paddingLeft) + parseFloat(style.paddingRight);
        var padY = parseFloat(style.paddingTop) + parseFloat(style.paddingBottom);
        mergedLayout.width = Math.max(rect.width - padX, 200);
        mergedLayout.height = Math.max(rect.height - padY, 300);

        Plotly.newPlot(elementId, plotData.traces || [], mergedLayout, {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d']
        });

        // Auto-resize on window resize
        if (!el._resizeAttached) {
            el._resizeAttached = true;
            window.addEventListener('resize', function() {
                if (el.querySelector('.js-plotly-plot')) {
                    Plotly.Plots.resize(el);
                }
            });
        }
    }

    function clearPlot(elementId) {
        var el = document.getElementById(elementId);
        if (!el) return;
        Plotly.purge(elementId);
    }

    window.Plots = {
        renderPlot: renderPlot,
        clearPlot: clearPlot,
        getDefaultLayout: getDefaultLayout
    };
})();
