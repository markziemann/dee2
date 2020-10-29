module KeyBoardHelpers exposing (..)
import Keyboard.Event exposing (KeyboardEvent)
import Maybe.Extra as MExtra

keyToMsg : String -> msg -> KeyboardEvent -> Maybe msg
keyToMsg keyCode msg keyboardEvent=
    case keyboardEvent.key of
        Just value ->
            if value == keyCode then
                Just msg
            else
                Nothing
        _ ->
            Nothing

try: List (KeyboardEvent -> Maybe msg) -> KeyboardEvent -> Maybe msg
try listeners keyboardEvent  =
        List.map (\func -> func keyboardEvent)
            listeners
        |> MExtra.orList

enterKey = keyToMsg "Enter"

arrowUp = keyToMsg "ArrowUp"

arrowDown = keyToMsg "ArrowDown"