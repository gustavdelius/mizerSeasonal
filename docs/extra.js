MathJax = {
  tex: {
    tags: 'ams',
    }
  };

  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml-full.js";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();
