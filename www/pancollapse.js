$(document).on('click', '.panel-heading-collapse', function(e){
    var panel = $(this).parents('.panel');
    var panelBody = panel.find('.panel-body');
	if(!panelBody.hasClass('panel-collapse')) {
		panelBody.slideUp();
        panelBody.addClass('panel-collapse collapse');
	} else {
        panelBody.slideDown();
        panelBody.removeClass('panel-collapse collapse');
	}
    panel.find('i').toggleClass('glyphicon-chevron-up glyphicon-chevron-down');
});