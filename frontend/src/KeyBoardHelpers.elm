module KeyBoardHelpers exposing (..)
import Keyboard.Event exposing (KeyboardEvent)

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

-- onKeyDown (considerKeyboardEvent ((arrowUp ArrowUp) |> orTy (arrowDown ArrowDown)))

orTry: (KeyboardEvent -> Maybe msg) -> (KeyboardEvent -> Maybe msg) -> KeyboardEvent -> Maybe msg
orTry key1 key2 keyBoardEvent=
    case key1 keyBoardEvent of
        Just value ->
            Just value
        Nothing ->
            case key2 keyBoardEvent of
                Just value ->
                    Just value
                Nothing ->
                    Nothing




enterKey = keyToMsg "Enter"

arrowUp = keyToMsg "ArrowUp"

arrowDown = keyToMsg "ArrowDown"