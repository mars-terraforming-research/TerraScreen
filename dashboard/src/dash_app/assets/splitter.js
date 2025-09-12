document.addEventListener('DOMContentLoaded', () => {
  const grid = document.getElementById('grid');
  if (!grid) return;
  const s1 = document.getElementById('split-1');
  const s2 = document.getElementById('split-2');
  let dragging = 0;

  function rect(){ return grid.getBoundingClientRect(); }
  function setLeft(px){ grid.style.setProperty('--left', Math.max(0, px) + 'px'); }
  function setRight(px){ grid.style.setProperty('--right', Math.max(0, px) + 'px'); }

  function onMove(clientX){
    if (!dragging) return;
    const r = rect();
    if (dragging === 1) setLeft(clientX - r.left);
    if (dragging === 2) setRight(r.right - clientX);
  }

  function onPointerMove(e){
    const x = (e.touches && e.touches[0]) ? e.touches[0].clientX : e.clientX;
    onMove(x);
    if (dragging) e.preventDefault();
  }

  function onPointerUp(){
    dragging = 0;
    document.body.style.cursor = '';
  }

  function start(which){ return (e) => {
    dragging = which;
    document.body.style.cursor = 'col-resize';
    e.preventDefault();
  };}

  s1.addEventListener('mousedown', start(1));
  s2.addEventListener('mousedown', start(2));
  s1.addEventListener('touchstart', start(1), {passive:false});
  s2.addEventListener('touchstart', start(2), {passive:false});

  window.addEventListener('mousemove', onPointerMove, {passive:false});
  window.addEventListener('touchmove', onPointerMove, {passive:false});
  window.addEventListener('mouseup', onPointerUp);
  window.addEventListener('touchend', onPointerUp);
});
