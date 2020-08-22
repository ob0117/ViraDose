
// dictionary for AAV Vector values
var VectorValues = {

    "AAV-1" : [4.47, 5.46],
    "AAV-2" : [2.00, 4.95],
    "AAV-4" : [3.81, 14.34],
    "AAV-5" : [3.65, 6.05],
    "AAV-6" : [6.65, 5.68],
    "AAV-7" : [3.60, 6.73],
    "AAV-8" : [14.55, 13.44],
    "AAV-9" : [2.23, 9.50]
};

// CWSF toggle button
function showCWSFDisplay() {

    var pdf = document.getElementById("pdf");

    if(pdf.style.display == "none") {
        pdf.style.display = "block";
    }
    else {
        pdf.style.display = "none";
    }
}

// Some values change based on AAV Vector
function changeAAVValues() {


    var AAVSelect = document.getElementById("aav-select");
    var vector = AAVSelect.options[AAVSelect.selectedIndex].value;

    var t = document.getElementById("t");
    var logc = document.getElementById("logc");

    t.value = VectorValues[vector][0];
    logc.value = VectorValues[vector][1];

}

// toggle math model select
function showMathModelFields() {

    var mathSelect = document.getElementById("math-select");
    var modelSelected = mathSelect.options[mathSelect.selectedIndex].value;

    var model1 = document.getElementById("model1params");
    var model2 = document.getElementById("model2params");

    if(modelSelected == "View innate immune response") {
        model1.style.display = "block";
        model2.style.display = "none";
    }
    else {
        model1.style.display = "none";
        model2.style.display = "block";
    }
}

// toggle advanced settings
function showAdvancedSettings() {

    advSettings = document.getElementById("adv-settings");

    if(advSettings.style.display == "none") {
        advSettings.style.display = "block";
    }
    else {
        advSettings.style.display = "none";
    }

}