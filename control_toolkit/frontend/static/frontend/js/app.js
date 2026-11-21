(function() {
    'use strict';

    var debounceTimer = null;

    function formatNumber(v) {
        if (v === null || v === undefined) return '—';
        if (typeof v === 'string') {
            var lower = v.toLowerCase();
            if (lower === 'nan' || lower === 'infinity' || lower === '-infinity' || lower === 'inf' || lower === '-inf') {
                return v;
            }
            v = parseFloat(v);
        }
        if (isNaN(v)) return 'NaN';
        if (!isFinite(v)) return v > 0 ? '∞' : '-∞';
        if (v === 0) return '0';
        return parseFloat(v.toPrecision(4)).toString();
    }

    function getSystem() {
        var numInput = document.getElementById('num-input');
        var denInput = document.getElementById('den-input');
        var unityCheckbox = document.getElementById('unity-checkbox');
        return {
            num: numInput ? numInput.value.trim() : '',
            den: denInput ? denInput.value.trim() : '',
            unity: unityCheckbox ? unityCheckbox.checked : false
        };
    }

    async function apiCall(endpoint, data) {
        var response = await fetch('/api/' + endpoint + '/', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(data)
        });
        if (!response.ok) {
            var errorData;
            try {
                errorData = await response.json();
            } catch (e) {
                throw new Error('Server error (' + response.status + ')');
            }
            throw new Error(errorData.error || errorData.detail || 'Server error (' + response.status + ')');
        }
        return await response.json();
    }

    function showToast(message, type) {
        type = type || 'error';
        var container = document.getElementById('toast-container');
        if (!container) return;
        var toast = document.createElement('div');
        toast.className = 'toast toast-' + type;
        toast.textContent = message;
        container.appendChild(toast);
        setTimeout(function() {
            toast.style.opacity = '0';
            setTimeout(function() {
                if (toast.parentNode) toast.parentNode.removeChild(toast);
            }, 300);
        }, 4000);
    }

    function showLoading(elementId) {
        var el = document.getElementById(elementId);
        if (!el) return;
        var spinner = document.createElement('div');
        spinner.className = 'loading-spinner';
        spinner.innerHTML = '<div class="spinner"></div>';
        el.appendChild(spinner);
    }

    function hideLoading(elementId) {
        var el = document.getElementById(elementId);
        if (!el) return;
        var spinners = el.querySelectorAll('.loading-spinner');
        spinners.forEach(function(s) { s.parentNode.removeChild(s); });
    }

    function updateTFDisplay() {
        if (debounceTimer) clearTimeout(debounceTimer);
        debounceTimer = setTimeout(async function() {
            var system = getSystem();
            if (!system.num && !system.den) return;
            var display = document.getElementById('tf-display');
            if (!display) return;
            try {
                var data = await apiCall('transfer-function', system);
                display.innerHTML = data.display || '';
            } catch (e) {
                display.innerHTML = '<span class="error">Error: ' + e.message + '</span>';
            }
        }, 300);
    }

    function initTabs() {
        var tabButtons = document.querySelectorAll('.tab-btn');
        tabButtons.forEach(function(btn) {
            btn.addEventListener('click', function() {
                tabButtons.forEach(function(b) { b.classList.remove('active'); });
                btn.classList.add('active');

                var target = btn.getAttribute('data-tab');
                var contents = document.querySelectorAll('.tab-content');
                contents.forEach(function(c) { c.classList.remove('active'); });
                var targetEl = document.getElementById('tab-' + target);
                if (targetEl) targetEl.classList.add('active');
            });
        });
    }

    document.addEventListener('DOMContentLoaded', function() {
        initTabs();

        var numInput = document.getElementById('num-input');
        var denInput = document.getElementById('den-input');
        var unityCheckbox = document.getElementById('unity-checkbox');

        if (numInput) numInput.addEventListener('input', updateTFDisplay);
        if (denInput) denInput.addEventListener('input', updateTFDisplay);
        if (unityCheckbox) unityCheckbox.addEventListener('change', updateTFDisplay);

        updateTFDisplay();
    });

    function formatComplex(c) {
        if (c === null || c === undefined) return '—';
        if (typeof c === 'string') return c;
        if (typeof c !== 'object') return formatNumber(c);
        var re = c.real !== undefined ? c.real : 0;
        var im = c.imag !== undefined ? c.imag : 0;
        if (Math.abs(im) < 1e-10) return formatNumber(re);
        if (Math.abs(re) < 1e-10) return formatNumber(im) + 'j';
        var sign = im >= 0 ? ' + ' : ' − ';
        return formatNumber(re) + sign + formatNumber(Math.abs(im)) + 'j';
    }

    function formatValue(v) {
        if (v === null || v === undefined) return '—';
        if (typeof v === 'string') return v;
        if (typeof v === 'object' && !Array.isArray(v)) {
            if ('real' in v && 'imag' in v) return formatComplex(v);
            if ('s' in v && 'K' in v) return 's = ' + formatNumber(v.s) + ', K = ' + formatNumber(v.K);
            if ('K' in v && 'real' in v && 'imag' in v) return 'K=' + formatNumber(v.K) + ' → ' + formatComplex(v);
            if ('K' in v && 'w' in v) return 'K = ' + formatNumber(v.K) + ', ω = ' + formatNumber(v.w) + ' rad/s';
            if ('lo' in v && 'hi' in v) return formatNumber(v.lo) + ' < K < ' + formatNumber(v.hi);
            if ('value' in v) {
                var s = v.value;
                if (v.wn !== undefined) s += '  ωₙ=' + formatNumber(v.wn);
                if (v.zeta !== undefined) s += '  ζ=' + formatNumber(v.zeta);
                return s;
            }
            return JSON.stringify(v);
        }
        return formatNumber(v);
    }

    window.App = {
        apiCall: apiCall,
        getSystem: getSystem,
        showToast: showToast,
        showLoading: showLoading,
        hideLoading: hideLoading,
        formatNumber: formatNumber,
        formatComplex: formatComplex,
        formatValue: formatValue
    };
})();
