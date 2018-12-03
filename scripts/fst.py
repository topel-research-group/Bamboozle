#!/usr/bin/env python3

#	Input from pipeline -> this script -> output Fst statistics table and graph
#
#	Copyright (C) 2018 Vilma Canfjorden. vilma.canfjorden@gmail.com
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <https://www.gnu.org/licenses/>.


import sys
import subprocess
import argparse
import fnmatch
import os
import glob
import csv

##################################################################################
parser = argparse.ArgumentParser(prog="fst.py")
parser.add_argument("-c", "--clean", action="store_true", help="Remove some files")
args = parser.parse_args()
##################################################################################
current_directory = os.getcwd()
name = os.path.basename(current_directory)
merged_vcf_pop1 = name + '_merged_pop1.vcf.gz'
merged_vcf_pop2 = name + '_merged_pop2.vcf.gz'
names1 = 'name_1_list.txt'#
names2 =  'name_2_list.txt'#
indv_txt_pop1 = name + '_indv_names_pop1.txt'
indv_txt_pop2 = name + '_indv_names_pop2.txt'
population_list = 'pop_list.txt'#
all_pop_merged = 'all_pop_merged.vcf.gz'
fst_out = 'pop1_pop2'
fst_out_in = '../Populations/pop1_pop2.weir.fst'# 
fst_out_flt = 'tmp.pop1_pop2_flt.table'
fst_out_flt_results = 'tmp.pop1_pop2_flt_results.table'
fst_out_flt2_results = 'tmp.pop1_pop2_flt2_results.table'
fst_results_sorted = 'pop1_pop2_flt_results_sorted.table'
fst_results_sorted_csv = 'pop1_pop2_flt_results_sorted.csv'
path_for_plot = 'Fst_stats/'#
add = '../'
##################################################################################

# Perform Fst-statistics on gziped vcf-files
def main():
	directories = current_directory + '/*/*/Bcftools/*.snpeff_annotated.vcf.gz'
	file_list = glob.glob(directories)
	for f in file_list:
		cmd1 = ['bcftools', 'index', '-c', '-f', f]
		process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
		while process1.wait() is None:
			pass
		process1.stdout.close()

	# Make directory for the merged vcf-files for population1 and population2
	population_directory = os.path.join(current_directory, r'Populations')
	if not os.path.exists(population_directory):
		os.makedirs(population_directory)

	# Making a list of vcf-files that will be input to bcftools merge and then merge population1
	directories2 = current_directory + '/*_1/*/Bcftools/*.snpeff_annotated.vcf.gz'
	name_list1 = glob.glob(directories2)
	myfile = open("name_1_list.txt","w")
	for n1 in name_list1:
		myfile.write("%s\n" % n1)

	myfile.close()
	cmd2 = ['bcftools', 'merge', '-l', add+names1, '-Oz', '-o', merged_vcf_pop1]   
	process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, cwd='Populations')
	while process2.wait() is None:
		pass
	process2.stdout.close()

	# Making a list of vcf-files that will be input to bcftools merge and then merge population2
	directories3 = current_directory + '/*_2/*/Bcftools/*.snpeff_annotated.vcf.gz'
	name_list2 = glob.glob(directories3)
	myfile2 = open("name_2_list.txt","w")
	for n2 in name_list2:
		myfile2.write("%s\n" % n2)

	myfile2.close()
	cmd3 = ['bcftools', 'merge', '-l', add+names2, '-Oz', '-o', merged_vcf_pop2]   
	process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, cwd='Populations')
	while process3.wait() is None:
		pass
	process3.stdout.close()

	# Making a txt file of the names of the individuals in the populations that is needed for vcftools --wei-fst-pop
	# and indexing the merged files for population1 and population2
	for file in os.listdir('Populations'):
		if fnmatch.fnmatch(file, '*_merged_pop1.vcf.gz'):
			cmd4 = ['bcftools', 'index', '-c', '-f', merged_vcf_pop1]
			process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, cwd='Populations')
			while process4.wait() is None:
				pass
			process4.stdout.close()

			cmd5 = ('bcftools query -l %s > %s') % (merged_vcf_pop1, indv_txt_pop1) 
			process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE, shell=True, cwd='Populations')
			while process5.wait() is None:
				pass
			process5.stdout.close()

		elif fnmatch.fnmatch(file, '*_merged_pop2.vcf.gz'):
			cmd6 = ['bcftools', 'index', '-c', '-f', merged_vcf_pop2]
			process6 = subprocess.Popen(cmd6, stdout=subprocess.PIPE, cwd='Populations')
			while process6.wait() is None:
				pass
			process6.stdout.close()

			cmd7 = ('bcftools query -l %s > %s') % (merged_vcf_pop2, indv_txt_pop2)
			process7 = subprocess.Popen(cmd7, stdout=subprocess.PIPE, shell=True, cwd='Populations')
			while process7.wait() is None:
				pass
			process7.stdout.close()

	
	# Making a list of vcf-files that will be input to bcftools merge and then merge population1 and population2
	# to a "all_merged" vcf file, which will be the input file to vcftools --weir-fst-pop 
	directories4 = current_directory + '/Populations/*_merged_*.vcf.gz'
	pop_list = glob.glob(directories4)
	myfile3 = open("pop_list.txt","w")
	for p in pop_list:
		myfile3.write("%s\n" % p)

	myfile3.close()
	cmd8 = ['bcftools', 'merge', '-l', add+population_list, '-Oz', '-o', all_pop_merged] 
	process8 = subprocess.Popen(cmd8, stdout=subprocess.PIPE, cwd='Populations')
	while process8.wait() is None:
		pass
	process8.stdout.close()

	# Making directory for Fst-results, input-files to highcharts
	fst_directory = os.path.join(current_directory, r'Fst_stats')
	if not os.path.exists(fst_directory):
		os.makedirs(fst_directory)

	# Fst_statistics 
	for file in os.listdir('Populations'):
		if fnmatch.fnmatch(file, 'all_pop_merged.vcf.gz'):
			cmd9 = ['vcftools', '--gzvcf', all_pop_merged, '--weir-fst-pop', indv_txt_pop1, '--weir-fst-pop', indv_txt_pop2, '--out', fst_out]
			process9 = subprocess.Popen(cmd9, stdout=subprocess.PIPE, cwd='Populations')
			while process9.wait() is None:
				pass
			process9.stdout.close()

	# Filtering the resulting files from vcftools and making a new directory called Fst_stats with the resulting files,  
	# the csv-file will be the input file to high charts 
	for file in os.listdir('Populations'):
		if fnmatch.fnmatch(file, '*.weir.fst'):
			cmd10 = ('cat %s | grep -v "nan" > %s') % (fst_out_in, fst_out_flt)
			process10 = subprocess.Popen(cmd10, stdout=subprocess.PIPE, shell=True, cwd='Fst_stats')
			while process10.wait() is None:
				pass
			process10.stdout.close()

	# Removing the results below zero
	for file in os.listdir('Fst_stats'):
		if fnmatch.fnmatch(file, '*flt.table'):
			cmd11 = ("awk '{if ($3 >0) print}' %s > %s") % (fst_out_flt, fst_out_flt_results)
			process11 = subprocess.Popen(cmd11, stdout=subprocess.PIPE, shell=True, cwd='Fst_stats')
			while process11.wait() is None:
				pass
			process11.stdout.close()

			# Rearrange columns (if needed)
			cmd12 = ('''awk '{print $1 "\\t" $2 "\\t" $3}' %s > %s''') % (fst_out_flt_results, fst_out_flt2_results)
			process12 = subprocess.Popen(cmd12, stdout=subprocess.PIPE, shell=True, cwd='Fst_stats')
			while process12.wait() is None:
				pass
			process12.stdout.close()

			# Sorting the POS column (needed for x-axis in highcharts)
			cmd13 = ("cat %s | sort -n > %s") % (fst_out_flt2_results, fst_results_sorted)
			process13 = subprocess.Popen(cmd13, stdout=subprocess.PIPE, shell=True, cwd='Fst_stats')
			while process13.wait() is None:
				pass
			process13.stdout.close()

			# Making a csv-file
			cmd14 = ('cat %s | tr "\\t" ","  > %s') % (fst_results_sorted, fst_results_sorted_csv)
			process14 = subprocess.Popen(cmd14, stdout=subprocess.PIPE, shell=True, cwd='Fst_stats')
			while process14.wait() is None:
				pass
			process14.stdout.close()

	# Making a plot of the Fst results using highcharts, the output is a html file
	for file in os.listdir('Fst_stats'):
		if fnmatch.fnmatch(file, 'pop1_pop2_flt_results_sorted.csv'):
			with open(path_for_plot + 'pop1_pop2_flt_results_sorted.csv') as infile: 
				csv_infile = csv.reader(infile)
				next(csv_infile)
				csv_list = [] 

				for line in csv_infile:
					csv_list.append('\\n\\"' +line[0]+ '\\"' + ";" + line[1] + ";" + line[2])  

			file=open(path_for_plot + 'Fst_results.html', 'w')

			string_prefix = '''"\\"CHROM\\";\\"POS\\";\\"WEIR_AND_COCKERHAM_FST\\"''' 
 
			html_part1= '''<div id="highcharts-aac96e77-0ad0-4715-a815-b427470cf979"></div><script>
					(function(){ var files = ["https://code.highcharts.com/stock/highstock.js","https://code.highcharts.com/highcharts-more.js","https://code.highcharts.com/highcharts-3d.js","https://code.highcharts.com/modules/data.js","https://code.highcharts.com/modules/exporting.js","http://code.highcharts.com/modules/funnel.js","http://code.highcharts.com/modules/solid-gauge.js"],loaded = 0; if (typeof window["HighchartsEditor"] === "undefined") {window.HighchartsEditor = {ondone: [cl],hasWrapped: false,hasLoaded: false};include(files[0]);} else {if (window.HighchartsEditor.hasLoaded) {cl();} else {window.HighchartsEditor.ondone.push(cl);}}function isScriptAlreadyIncluded(src){var scripts = document.getElementsByTagName("script");for (var i = 0; i < scripts.length; i++) {if (scripts[i].hasAttribute("src")) {if ((scripts[i].getAttribute("src") || "").indexOf(src) >= 0 || (scripts[i].getAttribute("src") === "http://code.highcharts.com/highcharts.js" && src === "https://code.highcharts.com/stock/highstock.js")) {return true;}}}return false;}function check() {if (loaded === files.length) {for (var i = 0; i < window.HighchartsEditor.ondone.length; i++) {try {window.HighchartsEditor.ondone[i]();} catch(e) {console.error(e);}}window.HighchartsEditor.hasLoaded = true;}}function include(script) {function next() {++loaded;if (loaded < files.length) {include(files[loaded]);}check();}if (isScriptAlreadyIncluded(script)) {return next();}var sc=document.createElement("script");sc.src = script;sc.type="text/javascript";sc.onload=function() { next(); };document.head.appendChild(sc);}function each(a, fn){if (typeof a.forEach !== "undefined"){a.forEach(fn);}else{for (var i = 0; i < a.length; i++){if (fn) {fn(a[i]);}}}}var inc = {},incl=[]; each(document.querySelectorAll("script"), function(t) {inc[t.src.substr(0, t.src.indexOf("?"))] = 1; }); function cl() {if(typeof window["Highcharts"] !== "undefined"){var options={"chart":{"type":"column","zoomType":"x"},"series[0]":{"type":"column"},"title":{"text":"Weir and Cockerham Fst"},"subtitle":{"text":"VCFtools - v0.1.13"},"series":[{"turboThreshold":0,"_colorIndex":0,"_symbolIndex":0,"type":"line","marker":{"enabled":false},"colorByPoint":false}],"data":{"csv":'''

			html_part2 = '''","googleSpreadsheetKey":false,"googleSpreadsheetWorksheet":false},"pane":{"background":[]},"responsive":{"rules":[]},"legend":{"enabled":true,"floating":false},"xAxis":[{}],"yAxis":[{"max":1}],"tooltip":{"shared":true},"mapNavigation":{"enableMouseWheelZoom":true,"enableButtons":true,"enabled":true},"plotOptions":{"series":{"animation":false}}};/*
// Sample of extending options:
Highcharts.merge(true, options, {
    chart: {
        backgroundColor: "#bada55"
    },
    plotOptions: {
        series: {
            cursor: "pointer",
            events: {
                click: function(event) {
                    alert(this.name + " clicked\n" +
                          "Alt: " + event.altKey + "\n" +
                          "Control: " + event.ctrlKey + "\n" +
                          "Shift: " + event.shiftKey + "\n");
                }
            }
        }
    }
});
*/new Highcharts.Chart("highcharts-aac96e77-0ad0-4715-a815-b427470cf979", options);}}})();
</script> '''
			csv_string = html_part1 + string_prefix + ''.join(csv_list) + html_part2

			file.write(csv_string)
			file.close()
				

	# Removing tmp-files
	if args.clean:
		for textfile in os.listdir('.'):
			if fnmatch.fnmatch(textfile, 'name_*_list.txt'):
				os.remove(textfile)
			elif fnmatch.fnmatch(textfile, 'pop_list.txt'):
				os.remove(textfile)
		
		for fstfile in os.listdir('Fst_stats'):
			if fnmatch.fnmatch(fstfile, 'tmp.*'):
				os.remove(current_directory + '/Fst_stats/' + fstfile)
		

if __name__ == "__main__":
	main()
