// shows the div with id id and hides the rest of the divs inside div
function show(element,id,div) {
	var parentDiv = document.getElementById(div);

	var children = parentDiv.children;
	for( var i=0; i<children.length; i++ ) {
		if( children[i].id === id ) {
			children[i].className = "shown";
		} else {
			children[i].className = "hidden";
		}
	}

	children = element.parentNode.children;
	for( var i=0; i<children.length; i++ ) {
		if( children[i].className.indexOf("selected") != -1 ) {
			children[i].className = '';
		} 
	}
	element.className = "selected";
}