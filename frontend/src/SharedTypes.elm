module SharedTypes exposing (..)

import Http


type alias WebData a =
    RemoteData Http.Error a


unwrapWebData : b -> (a -> b) -> WebData a -> b
unwrapWebData default func webData =
    case webData of
        Success a ->
            func a

        _ ->
            default


toWebData : Result Http.Error a -> WebData a
toWebData result =
    case result of
        Ok value ->
            Success value

        Err err ->
            Failure err


type RemoteData e a
    = NotAsked
    | Loading
    | Failure e
    | Success a


type alias PaginationOffset =
    { perPage : Int
    , offset : Int
    }
