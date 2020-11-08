import {Elm} from './src/Main.elm'

var app = Elm.Main.init({
    node: document.querySelector('main')
});


app.ports.isElementTextTruncated.subscribe(function (id_element) {
    let id = id_element[0];
    let element = id_element[1];
    console.log(element.offsetWidth, element.scrollWidth);
    app.ports.receiveIdElement.send(
        [id, element.offsetWidth < element.scrollWidth]
    );
});
