module KeyBoardHelpers exposing (..)
import Keyboard.Event exposing (KeyboardEvent)
import Maybe.Extra as MExtra
import Keyboard.Key as KKey

keyToMsg : (String, KKey.Key) -> msg -> KeyboardEvent -> Maybe msg
keyToMsg (key, keyCode) msg keyboardEvent=
    let
        eventKey = Maybe.withDefault "" keyboardEvent.key
        eventKeyCode = keyboardEvent.keyCode
    in
    if (key == eventKey) || (keyCode == eventKeyCode) then
        Just msg
    else
        Nothing

try: List (KeyboardEvent -> Maybe msg) -> KeyboardEvent -> Maybe msg
try listeners keyboardEvent  =
        List.map (\func -> func keyboardEvent)
            listeners
        |> MExtra.orList

enterKey = keyToMsg ("Enter", KKey.Enter)

arrowUp = keyToMsg ("ArrowUp", KKey.Up)

arrowDown = keyToMsg ("ArrowDown", KKey.Down)