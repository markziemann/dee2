import {Elm} from './src/Main.elm'

var app = Elm.Main.init({
    node: document.querySelector('main')
});

app.ports.consoleLog.subscribe(function (message) {
    console.log(message)
});