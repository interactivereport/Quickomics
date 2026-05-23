(function () {
  var style = document.createElement('style');
  style.textContent = `
    .ui-sortable-handle.multidrag-selected {
      background-color: #d4edda !important;
      border: 2px solid #28a745 !important;
    }
    .ui-sortable-placeholder {
      display: none !important;
      height: 0px !important;
      visibility: hidden !important;
    }
  `;
  document.head.appendChild(style);
})();

$(document).on('shiny:connected', function () {

  // Ctrl/Cmd+click to toggle selection
  $(document).on('click', '.ui-sortable-handle', function (e) {
    if (e.ctrlKey || e.metaKey) {
      e.preventDefault();
      $(this).toggleClass('multidrag-selected');
    } else {
      $('.ui-sortable-handle.multidrag-selected').not(this).removeClass('multidrag-selected');
    }
  });

  // Capture and hide companions when drag starts
  $(document).on('sortstart', function (event, ui) {
    var $list = $(event.target);

    if (!ui.item.hasClass('multidrag-selected')) {
      $list.find('.multidrag-selected').removeClass('multidrag-selected');
    }

    var $companions = $list.find('.multidrag-selected').not(ui.item);
    ui.item.data('multidrag-companions', $companions);
    $companions.hide();
  });

  // Snap companions in place, then update Shiny for both lists
  $(document).on('sortstop', function (event, ui) {
    var $companions = ui.item.data('multidrag-companions') || $();
    if ($companions.length) {
      $companions.insertAfter(ui.item).show();
    }
    $('.multidrag-selected').removeClass('multidrag-selected');

    // Update this list and its connected partner
    var $list = $(event.target);
    var connectedSelector = $list.sortable('option', 'connectWith');

    [$list, $(connectedSelector)].forEach(function ($l) {
      if (!$l || !$l.length) return;
      var id = $l.attr('id');
      if (id) {
        Shiny.setInputValue(
          id,
          $l.sortable('toArray', { attribute: 'data-value' }),
          { priority: 'event' }
        );
      }
    });
  });

});