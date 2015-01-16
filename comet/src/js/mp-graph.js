function mpgraph(container, nodes, links){
	// Style and formatting
	var margin = {left: 10, right: 40, top: 10, bottom: 10},
		width = 315,
		height = 120 + 20 * (Math.max(0, nodes.length-6)),
		format = d3.format(".2g");

	var color = d3.scale.category20();

	// Set up the force layout
	var force = d3.layout.force()
		.charge(-120)
		.linkDistance(30)
		.size([width, height])
		.on("tick", tick);

	var drag = force.drag()
		.on("dragstart", dragstart)
		.on("dragend", function(){ force.stop(); });

	// Make sure the nodes are fixed to begin with, and update the locations
	// to match the size of the plot
	nodes.forEach(function(d){
		d.x = margin.left + d.x * (width - margin.right - margin.left);
		d.y = margin.top  + d.y * (height - margin.top - margin.bottom);
		d.fixed  = true; // means the force simulation won't move them around
	});

	// Initalize the force layout
	force.nodes(nodes)
		.links(links)
		.start()
		.stop();

	// Add the SVG
	var svg = container.append("svg")
		.attr("width", width)
		.attr("height", height);

	// Add links in groups, so edge labels stay with lines
	var link = svg.selectAll(".link")
		.data(links).enter()
		.append("line")
		.attr("class", "link")
		.style("stroke", "#000")
		.style("stroke-width", function(d) { return Math.sqrt(d.weight); });

	var edgeLabels = svg.selectAll("text")
		.data(links).enter()
		.append("text")
		.style("font-size", 8)
		.text(function(d) { return format(d.weight); });

	// Add nodes in groups, so that the label stays with the node
	var node = svg.selectAll(".node")
		.data(nodes).enter()
		.append("g")
		.call(drag);

	node.append("circle")
		.attr("class", "node")
		.attr("r", 10)
		.style("fill", function(){ return color(Math.random() * 20); })

	node.append("title")
		.text(function(d) { return d.name; });

	node.append("text")
		.attr({"dx":10,"dy":5,"class":"nodelabel","font-size": 12 })
		.text(function(d){return d.name;});

	function tick() {
		// Move the edges, edgeLabels, and nodes into place
		link.attr("x1", function(d) { return d.source.x; })
			.attr("y1", function(d) { return d.source.y; })
			.attr("x2", function(d) { return d.target.x; })
			.attr("y2", function(d) { return d.target.y; });

		edgeLabels.attr("x", function(d) { return d.source.x + (d.target.x - d.source.x) / 2; })
			.attr("y", function(d) { return d.source.y + (d.target.y - d.source.y) / 2; });

		node.attr("transform", function(d){ return "translate(" + d.x + "," + d.y + ")"; })
	}

	function dragstart(d) {
		d3.select(this).select("circle").classed("fixed", d.fixed = true);
	}

	tick(); // get everything started
}