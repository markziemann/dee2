module Helpers exposing (..)

import Array
import Http exposing (Error(..))
import Json.Decode as Decode exposing (Decoder, field)
import SearchPage.Types exposing (SearchResults)
import SharedTypes


decodeSearchResults : SharedTypes.PaginationOffset -> Decoder SearchResults
decodeSearchResults ({ offset } as paginationOffset) =
    Decode.map2 SearchResults
        (field "hits" Decode.int)
        (field "rows"
            (Decode.list (Decode.dict Decode.string)
                -- Decode.map doesn't iterate (confusing!) its more like function application
                -- it should be called apply
                |> Decode.map Array.fromList
                |> Decode.map (Array.indexedMap (\idx data -> SearchPage.Types.SearchResult (idx + offset) data))
            )
        )

errorToString : Http.Error -> String
errorToString error =
    case error of
        BadUrl url ->
            "The URL " ++ url ++ " was invalid"

        Timeout ->
            "Unable to reach the server, try again"

        NetworkError ->
            "Unable to reach the server, check your network connection"

        BadStatus 500 ->
            "The server had a problem, try again later"

        BadStatus 400 ->
            "Verify your information and try again"

        BadStatus 404 ->
            "Failed to find resource"

        BadStatus 503 ->
            "Service Unavailable"

        BadStatus _ ->
            "Unknown Error"

        BadBody errorMessage ->
            errorMessage
