// Function to oggle the display style of an element
function toggleBlock(id,button_id,show_text,hide_text) {
  var x = document.getElementById(id);
  var b = document.getElementById(button_id);
  if (x.style.display === "none") {
    x.style.display = "block";
    b.innerHTML = hide_text;
  } else {
    x.style.display = "none";
    b.innerHTML = show_text;
  }
}
