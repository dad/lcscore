function show(evt, id, subelements) {
    "use strict";
    for (var ei = 0; ei < subelements.length; ei++) {
    	var subelement = subelements[ei]
		var fullid = id.concat("-").concat(subelement);
		console.log(fullid);
		var obj = evt.target.ownerDocument.getElementById(fullid);
		//obj.visibility = 'visible'
		obj.setAttribute("display", "inline");
    }
}
function hide(evt, id, subelements) {
    "use strict";
    for (var ei = 0; ei < subelements.length; ei++) {
    	var subelement = subelements[ei]
		var fullid = id.concat("-").concat(subelement);
		console.log(fullid);
		var obj = evt.target.ownerDocument.getElementById(fullid);
		var tog = evt.target.getAttribute('show-toggle')
		if (tog == 0 || tog == null) {
			obj.setAttribute("display", "none");
		}
    }
}
function borderHighlight(evt) {
    "use strict";
	var tog = evt.target.getAttribute('show-toggle')
    if (tog == 0 || tog == null) {
	    evt.target.setAttribute('stroke-width','2.0');
	    evt.target.setAttribute('stroke','black');
	}
}
function borderRestore(evt) {
    "use strict";
	var tog = evt.target.getAttribute('show-toggle')
	if (tog == 0 || tog == null) {
		evt.target.setAttribute('stroke-width','1.0');
		evt.target.setAttribute('stroke','none');
	}
}
function toggleVisibility(evt) {
    "use strict";
    var val = evt.target.getAttribute('show-toggle')
    var newval;
    if (val == 1) {
    	evt.target.setAttribute('show-toggle',0)
   		evt.target.setAttribute('stroke-width','2.0');
    } else if (val == 0 || val == null) {
    	evt.target.setAttribute('show-toggle',1)
   		evt.target.setAttribute('stroke-width','3.0');
    }
	
}

