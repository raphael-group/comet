function dplusViz(data){
	// Extract the weights
	var weights = Object.keys(data.graphs).sort(d3.ascending).map(function(n){ return +n; });

	// Select the various elements on the page
	var slider = $("input#min-weight"),
		currentWeight = d3.select("span#weight-threshold"),
		results = d3.select("div#results");

	// Initialize the slider
	var initialWeight = 1.*location.hash.split("#")[1],
		initialIndex =  weights.indexOf(initialWeight);

	initialIndex = initialIndex == - 1 ? 1 : initialIndex;
	slider.attr("max", weights.length-1);
	slider.val(initialIndex);
	currentWeight.text(weights[initialIndex]);

	// Watch the slider to update on change
	slider.on("change", drawComponents);

	// Add the statistical significance (if applicable)
	if (data.stats){
		// Create a list of the delta values (not including p-values) for
		// easy lookup later
		data.stats.deltaVals = data.stats.deltas.map(function(d){ return d.delta; });

		// Initialize the slider to start at the first significant result
		slider.val(weights.indexOf(data.stats.deltaVals[0]));
		currentWeight.text(data.stats.deltaVals[0]);

		// Add the significance plot
		d3.select("h4#significance-title").text("Significance Plot (N=" + data.stats.N + ")")
		d3.select("#significance").style("display", "inline");
		d3.select("#signifiance-plot").html(data.stats.plot);

		// Update the select with different options for deltas that have
		// statistical significance
		var deltaSelect = d3.select("#choose-stat-delta select"),
			deltaSelectJQ = $("#choose-stat-delta select"); // jquery version so we can access the value

		d3.select("#choose-stat-delta").style("display", "inline");
		deltaSelect.selectAll(".opt")
			.data(data.stats.deltas).enter()
			.append("option")
			.attr("value", function(d){ return d.delta; })
			.text(function(d){ return d.delta +  " (P=" + d.pval + ")"; })

		// Update the slider with the value from the select and draw the
		// updated components
		deltaSelect.on("change", function(){
			slider.val(weights.indexOf(+deltaSelectJQ.val()));
			drawComponents();
		})
	}

	// Set up the crossfilter
	var collections = crossfilter(data.tables);
	var geneDimension = collections.dimension(function(d){ return d.allGenes; })

	function drawComponents(){
		// Extract the current minimum edge weight
		var index      = slider.val(),
			weight     = weights[index],
			components = data.graphs[weight];

		// Update the slider and statistical significance select (if applicable)
		currentWeight.text(weight);
		location.hash = weight;
		if (data.stats){
			var val = data.stats.deltaVals.indexOf(weight) == -1 ? "-1" : weight;
			deltaSelectJQ.val(val);
		}

		// Add the marginal probability graph, mutation matricx button, and
		// sampled gene set button
		results.selectAll("div.cc").remove();
		components.forEach(function(obj, index){
			// Copy the links/nodes so they can be modified
			var links = obj.edges.map(function(d){ return $.extend(true, {}, d); }),
				nodes = obj.nodes.map(function(d){ return $.extend(true, {}, d); }),
				nodeNames = nodes.map(function(d){ return d.name; });

			// Add a container to hold the graphs
			var container = results.append("div")
				.attr("class", "cc list-group-item col-lg-12")

			// Add a column of graphs
			var graphContainer = container.append("div")
				.attr("class", "col-lg-4")
			mpgraph(graphContainer, nodes, links);

			// Add a column of mutation matrices
			var mtxContainer = container.append("div")
				.attr("class", "col-lg-8");

			mtxContainer.append("button")
				.text("Draw mutation matrix")
				.style("margin", "0px auto")
				.on("click", function(){
				    d3.select(this).remove();
				    var sampleToTypes = filterSamples(nodes, data.mutations.M, data.mutations.sampleToTypes),
				    	mutationData  = {M :{}, sampleToTypes: sampleToTypes, typeToSamples: data.mutations.typeToSamples};
				    nodes.forEach(function(d){
						mutationData.M[d.name] = data.mutations.M[d.name];
				    });
				    var m2Style = { style: { width: 690, labelHeight: 0 } };
					var m2Chart = mutation_matrix(m2Style).addCoverage();

					mtxContainer.datum(mutationData);
					m2Chart(mtxContainer);
			});

			// Add a new row for the table of results
			var tbl = container.append("div").attr("class", "col-lg-12");

			tbl.selectAll(".btn")
				.data([{visible: false}]).enter()
				.append("button")
				.text("Show sampled collections of gene sets")
				.style("margin-bottom", "10px")
				.on("click", function(d){
					if (d.visible){
						d3.select(this).text("Show sampled collections of gene sets");
						hideTable();
					}
					else{
						d3.select(this).text("Hide sampled collections of gene sets");
						showTable();
					}
					d.visible = d.visible ? false : true;
				});

			function hideTable(){
				tbl.select("table").remove();
				tbl.selectAll("div").remove();
			}

			function showTable(d){
					// Extract relevant rows
					var tableRows = geneDimension.filter(function(d){
						return nodes.reduce(function(s, n){
							return d[n.name] ? s + 1 : s;
						}, 0) > 1;
					}).top(Infinity);

					// Add a table
					var T = tbl.append("table")
						.attr("id", "table" + index)
						.attr("class", "table table-striped")
						.style("font-size", "10px");

					// Create a list of column headers
					var headerNames = [{title: "Freq", orderable: true}, {title: "W", orderable: true}];
					headerNames = headerNames.concat(d3.range(tableRows[0].genesets.length*2).map(function(i){
						return i % 2 ? {title: "&phi;(M)", orderable: true} : {title: "Genes", orderable: false};
					}));

					// Create a 2D array of rows
					var rows = [];
					tableRows.forEach(function(d){
						var row = [d.freq, d.totalWeight];
						d.genesets.forEach(function(G){
							row.push(G.genes.map(function(g){
								// Highlight the nodes in the gene set in yellow
								return nodeNames.indexOf(g) == -1 ? g : "<span style='background:yellow'>" + g + "</span>";
							}).join(", "));
							row.push(G.pval);
						});
						rows.push(row);
					});

					// Add the dataTable
					$("#table" + index).dataTable({
						"data": rows,
						"columns": headerNames.map(function(d){ return {title: d.title, orderable: d.orderable}; }),
						"order": [[0, "desc"], [1, "desc"]],
						"lengthChange": false,
						"deferRender": true,
						"searching": false
					});
				}
		});	
	}
	// end drawComponents
	drawComponents();

}

function filterSamples(genes, M, sampleToTypes){
    // Make a marking samples that are mutated
    var sampleToMutated = {};
    genes.forEach(function(g){
        Object.keys(M[g.name]).forEach(function(s){
            sampleToMutated[s] = true;
        });
    });

    // Create a new sample-to-type dictionary containing only
    // those samples mutated in the gene set
    var filteredSampleToTypes = {};
    Object.keys(sampleToMutated).forEach(function(s){
        filteredSampleToTypes[s] = sampleToTypes[s];
    });

    return filteredSampleToTypes;
}