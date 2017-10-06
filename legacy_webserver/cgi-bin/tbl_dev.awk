awk 'BEGIN { print "<table border="1"><tr><th> <input type=\"checkbox\" name=\"DataSetList\" onClick=\"toggle(this)\" />Select all</th><th> SRA run accession</th><th> QC summary </th><th>SRA experiment accession</th><th>SRA sample accession</th><th>SRA project accession</th><th>SRA submission accession</th><th>GEO series accession</th><th>GEO sample accession</th></tr>" }
	{ print "<tr><td> <input type='checkbox' name='x' value="$1">  </td><td>  <a href=http://www.ncbi.nlm.nih.gov/sra/?term="$1" >"$1"</a>  </td><td>  <a href=/qc/"$1".qc > "$2" </a> </td><td>" $3 "</td><td>" $4 "</td><td>" $5 "</td><td>" $6 "</td><td>" $7 "</td><td>" $8 "</td></tr>" }
     END   { print "</table>" }'

