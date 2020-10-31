module KeyBoardHelpers exposing (..)

import Json.Decode exposing (..)


arrowUp =
    maybeKeyEvent "ArrowUp"


arrowDown =
    maybeKeyEvent "ArrowDown"


enter =
    maybeKeyEvent "Enter"


maybeKeyEvent : String -> msg -> Decoder msg
maybeKeyEvent str msg =
    andThen
        (\result ->
            case result == str of
                True ->
                    succeed msg

                False ->
                    fail "Ignoring keyboard event"
        )
        (field "key" string)


