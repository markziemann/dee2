/*******************************************************************************
 SWIM2.0 :: Simple website menu
 ------------------------------------------------------------------------------
 Copyright (c) 2005 James Edwards (brothercake)          <cake@brothercake.com>
 BSD License                          See license.txt for licensing information
 Info/Docs                         http://www.brothercake.com/scripts/listmenu/
 ------------------------------------------------------------------------------
 For professional menu solutions visit                     http://www.udm4.com/ 
 ------------------------------------------------------------------------------
*******************************************************************************/



window.onload = function()
{
	var verticals = new simpleMenu('menu-v', 'vertical');
	var horizontals = new simpleMenu('menu-h', 'horizontal');
};



function simpleMenu(navid, orient){if(typeof document.getElementById == 'undefined' || /opera[\/ ][56]/i.test(navigator.userAgent)) { return; }this.iskde = navigator.vendor == 'KDE';this.isie = typeof document.all != 'undefined' && typeof window.opera == 'undefined' && !this.iskde;this.isoldsaf = navigator.vendor == 'Apple Computer, Inc.' && typeof XMLHttpRequest == 'undefined';this.tree = document.getElementById(navid);if(this.tree != null){this.items = this.tree.getElementsByTagName('li');this.itemsLen = this.items.length;var i = 0; do{this.init(this.items[i], this.isie, this.isoldsaf, this.iskde, navid, orient);}while (++i < this.itemsLen);}}simpleMenu.prototype.init = function(trigger, isie, isoldsaf, iskde, navid, ishoriz){trigger.menu = trigger.getElementsByTagName('ul').length > 0 ? trigger.getElementsByTagName('ul')[0] : null;trigger.link = trigger.getElementsByTagName('a')[0];trigger.issub = trigger.parentNode.id == navid;trigger.ishoriz = ishoriz == 'horizontal';this.openers = { 'm' : 'onmouseover', 'k' : (isie ? 'onactivate' : 'onfocus') };for(var i in this.openers){trigger[this.openers[i]] = function(e){if(!iskde) { trigger.link.className += (trigger.link.className == '' ? '' : ' ') + 'rollover'; }if(trigger.menu != null){if(trigger.ishoriz) { trigger.menu.style.left = (isie || isoldsaf) ? trigger.offsetLeft + 'px' : 'auto'; }trigger.menu.style.top = (trigger.ishoriz && trigger.issub) ? (isie || (trigger.ishoriz && isoldsaf)) ? trigger.link.offsetHeight + 'px' : 'auto' : (isie || (trigger.ishoriz && isoldsaf)) ? trigger.offsetTop + 'px' : '0';}};}this.closers = { 'm' : 'onmouseout', 'k' : (isie ? 'ondeactivate' : 'onblur') };for(i in this.closers){trigger[this.closers[i]] = function(e){this.related = (!e) ? window.event.toElement : e.relatedTarget;if(!this.contains(this.related)){if(!iskde) { trigger.link.className = trigger.link.className.replace(/[ ]?rollover/g, ''); }if(trigger.menu != null){trigger.menu.style[(trigger.ishoriz ? 'left' : 'top')] = trigger.ishoriz ? '-10000px' : '-100em';}}};}if(!isie){trigger.contains = function(node){if (node == null) { return false; }if (node == this) { return true; }else { return this.contains(node.parentNode); }};}}
