var width = 760,
    height = 500

var svg = d3.select("#pics").append("svg")
    .attr("width", width)
    .attr("height", height);

var force = d3.layout.force()
    .gravity(.15)
    .distance(30)
    .charge(-80)
    .size([width, height]);

var dat = d3.json("../data/graph.json", function(json){dat=json});
d3.json("../data/members.json", function(json) {
  force
      .nodes(json.results)
      .start();

  var link = svg.selectAll(".link")
      .data(dat.links)
    .enter().append("line")
      .attr("class", "link");

  var node = svg.selectAll(".node")
      .data(json.results)
    .enter().append("g")
      .attr("class", "node")
      .call(force.drag);

  node.append("image")
      .attr("xlink:href",function(d){return d.photo_url})
      .attr("x", -8)
      .attr("y", -8)
      .attr("width", 50)
      .attr("height", 50)
		.on("mouseover", function() {
        d3.select(this).transition()
        .attr('width', 150)
        .attr('height',150)
        .delay(0)
        .duration(500)
		.ease("elastic", 10, 1)         
    })
      .on("mouseout", function() {
        d3.select(this).transition()
        .attr('width', 50)
        .attr('height',50)
        .delay(0)
        .duration(500)});

  node.append("text")
      .attr("dx", 12)
      .attr("dy", ".35em")
      .text(function(d) { return d.name });

  force.on("tick", function() {
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });
  });
});
