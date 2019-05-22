#!/usr/bin/env python3

# This script makes a highcharts plot of the file provided.
# The input file should be a csv file with headers "CHROM,POS,VALUE"
# If the file is big use the `make_plot.py` utility script instead 
# by using the same input file but only plotting a selected contig.  
# Or big files can be plotted using the `pandas_fst_plot.py`
# (a non-interactive plot).

import csv
import argparse
import subprocess
import os
import fnmatch

#######################################################################

parser = argparse.ArgumentParser(prog="make_plot.py")
parser.add_argument("-i", "--inputfile", \
		help="input file, must be csv file", required=True)
parser.add_argument("-y", "--yaxis", \
		help="Name of header for yaxis", \
		required=True)
args = parser.parse_args()

#######################################################################

filtered = 'tmp_filtered.csv'
infile = args.inputfile 
yaxis = args.yaxis
name = os.path.basename(infile)
output_name = name + "_results.html"

#######################################################################

def main():
	with open(infile) as i:
		csv_filtered = csv.reader(i)
		next(csv_filtered, None)
		csv_list = []
		for line in csv_filtered:
			csv_list.append('\\n\\"' + line[0] + '\\"' + ";" + line[1] + ";" + line[2])

	file=open(output_name, 'w')
	html_part1= '''<div id="highcharts-aac96e77-0ad0-4715-a815-b427470cf979"></div><script>
(function(){ var files = ["https://code.highcharts.com/stock/highstock.js","https://code.highcharts.com/highcharts-more.js","https://code.highcharts.com/highcharts-3d.js","https://code.highcharts.com/modules/data.js","https://code.highcharts.com/modules/exporting.js","http://code.highcharts.com/modules/funnel.js","http://code.highcharts.com/modules/solid-gauge.js"],loaded = 0; if (typeof window["HighchartsEditor"] === "undefined") {window.HighchartsEditor = {ondone: [cl],hasWrapped: false,hasLoaded: false};include(files[0]);} else {if (window.HighchartsEditor.hasLoaded) {cl();} else {window.HighchartsEditor.ondone.push(cl);}}function isScriptAlreadyIncluded(src){var scripts = document.getElementsByTagName("script");for (var i = 0; i < scripts.length; i++) {if (scripts[i].hasAttribute("src")) {if ((scripts[i].getAttribute("src") || "").indexOf(src) >= 0 || (scripts[i].getAttribute("src") === "http://code.highcharts.com/highcharts.js" && src === "https://code.highcharts.com/stock/highstock.js")) {return true;}}}return false;}function check() {if (loaded === files.length) {for (var i = 0; i < window.HighchartsEditor.ondone.length; i++) {try {window.HighchartsEditor.ondone[i]();} catch(e) {console.error(e);}}window.HighchartsEditor.hasLoaded = true;}}function include(script) {function next() {++loaded;if (loaded < files.length) {include(files[loaded]);}check();}if (isScriptAlreadyIncluded(script)) {return next();}var sc=document.createElement("script");sc.src = script;sc.type="text/javascript";sc.onload=function() { next(); };document.head.appendChild(sc);}function each(a, fn){if (typeof a.forEach !== "undefined"){a.forEach(fn);}else{for (var i = 0; i < a.length; i++){if (fn) {fn(a[i]);}}}}var inc = {},incl=[]; each(document.querySelectorAll("script"), function(t) {inc[t.src.substr(0, t.src.indexOf("?"))] = 1; }); function cl() {if(typeof window["Highcharts"] !== "undefined"){var options={"chart":{"type":"line","zoomType":"x"},"series[0]":{"type":"line"},"title":{"text":""},"subtitle":{"text":"input: %s"},"series":[{"turboThreshold":0,"_colorIndex":0,"_symbolIndex":0,"type":"line","marker":{"enabled":false},"colorByPoint":false}],"data":{"csv":'''%(name)

	string_prefix = '''"\\"CHROM\\";\\"POS\\";\\"%s\\"'''%(yaxis)

	html_part2 = '''","googleSpreadsheetKey":false,"googleSpreadsheetWorksheet":false},"pane":{"background":[]},"responsive":{"rules":[]},"legend":{"enabled":true,"floating":false},"xAxis":[{}],"yAxis":[{"max":1}],"tooltip":{"shared":true},"mapNavigation":{"enableMouseWheelZoom":true,"enableButtons":true,"enabled":true},"plotOptions":{"series":{"animation":false}}};new Highcharts.Chart("highcharts-aac96e77-0ad0-4715-a815-b427470cf979", options);}}})();
</script> '''
	csv_string = html_part1 + string_prefix + ''.join(csv_list) + html_part2
	file.write(csv_string)
	file.close()

	for fstfile in os.listdir('.'):
		if fnmatch.fnmatch(fstfile, 'tmp_filtered.csv'):
			os.remove(fstfile)

if __name__ == "__main__":
	main()
