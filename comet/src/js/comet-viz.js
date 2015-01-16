function dplusViz(data){
	// Extract the weights
	var graph = data.graph,
		numEdges = graph.edges.length,
		weights = graph.weights;

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
		var index  = slider.val(),
			weight = weights[index];

		// Update the slider and statistical significance select (if applicable)
		currentWeight.text(weight);
		location.hash = weight;
		if (data.stats){
			var val = data.stats.deltaVals.indexOf(weight) == -1 ? "-1" : weight;
			deltaSelectJQ.val(val);
		}

		// Add/Remove the edges as needed, then extract the current components
		var edges = graph.edges.filter(function(d){ return d.weight >= weight; }),
			components = connectedComponents(edges, 2);

		// Add the marginal probability graph, mutation matricx button, and
		// sampled gene set button
		results.selectAll("div.cc").remove();
		components.forEach(function(cc, index){
			// Copy the links/nodes so they can be modified
			var links = edges.filter(function(d){ return cc.indexOf(d.source) !== -1; })
							 .map(function(d){ return $.extend(true, {}, d); }),
				nodes = cc.map(function(i){ return $.extend(true, {}, graph.nodes[i]); }),
				nodeNames = nodes.map(function(d){ return d.name; });

			// Renumber the edges from 0..len(component)
			var nodeToIndex = {};
			cc.forEach(function(n, i){ nodeToIndex[n] = i; });
			links.forEach(function(d){
				d.source = nodeToIndex[d.source];
				d.target = nodeToIndex[d.target];
			});

			// Compute a position for each node in a circular layout
			// (based off of http://goo.gl/gqvhgt).
			var increase = Math.PI * 2 / nodes.length,
				x = 0,
				y = 0,
				angle = 0;

			nodes.forEach(function(n, i){
				n.x = 0.5 * Math.cos(angle) + 0.5;
				n.y = 0.5 * Math.sin(angle) + 0.5;
				angle += increase;
			});

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

function connectedComponents(edges, minCCSize){
	UF = UnionFind();
	edges.forEach(function(d){ UF.union([d.source, d.target]); });
	return UF.groups().filter(function(g){ return g.length >= minCCSize; });
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

// Implementation of UnionFind data structure
// Based off of NetworkX's Python implementation: https://networkx.github.io/documentation/latest/_modules/networkx/utils/union_find.html
function UnionFind(){
	// Instance variables
	var weights = {},
		parents = {};

	// Find and return the name of the set containing the object
	function get(x){
		// check for previously unknown object
		if (!(x in parents)){
			parents[x] = x;
			weights[x] = 1;
			return x;
		} else {
			// find path of objects leading to the root
			var path = [x],
				root = parents[x],
				count = 0;

			while (root != path[path.length - 1] && count <= 15){
				path.push( root );
				root = parents[root];
				count++;
			}

			// compress the path and return
			path.forEach(function(ancestor){
				parents[ancestor] = root;
			});

			return root;
		}
	}

	// Find the sets containing the objects and merge them all
	function union(xs){
		// Convert xs to a list if it isn't one already
		if (xs.constructor != Array){
			xs = [xs];
		}

		// Merge all sets containing any x in xs
		var roots = xs.map(get),
			heaviest = d3.max(roots.map(function(r){ return [weights[r], r]; }))[1];

		roots.forEach(function(r){
			if (r != heaviest){
				weights[heaviest] += weights[r];
				parents[r] = heaviest;
			}
		});
	}

	// Return a list of lists containing each group
	function groups(){
		var groupIndex = 0,
			groupToIndex = {},
			currentGroups = [[]];

		Object.keys(parents).forEach(function(n){
			var group = get(n);
			if (!(group in groupToIndex)) groupToIndex[group] = groupIndex++;
			if (currentGroups.length <= groupToIndex[group]) currentGroups.push([]);
			currentGroups[groupToIndex[group]].push( +n );
		});

		return currentGroups;
	}

	return { get: get, union: union, groups: groups };
}