module MainViews exposing (..)

import Array
import Dict
import Html exposing (..)
import Html.Attributes as Attr
import Html.Events as Events
import MainTypes exposing (..)
import SearchBarTypes exposing (SearchResult, SearchResults)
import Table


selectClickedResult : SearchResult -> List (Html.Attribute Msg)
selectClickedResult result =
    [ Attr.class
        (if result.selected then
            "table-primary"

         else
            ""
        )

    --    Probably shouldn't be passing functions to update
    , ResultClicked result
        |> Events.onClick
    ]


listWrapped a =
    (::) a []



--viewSearchResults : SearchResults -> Html Msg
--viewSearchResults searchResults =
--    searchResults
--        |> Array.map
--            (\result ->
--                tr (selectClickedResult result)
--                    (List.map (\( key, value ) -> td [] [ text value ]) result.data)
--            )
--        |> Array.toList
--        |> tbody []
--        |> listWrapped
--        |> table [ Attr.class "table table-hover table-sm table-bordered table-responsive" ]
--type alias SearchResult =
--    { id : Int
--    , data : SearchData -- Array (Dict String String)
--    , selected : Bool
--    }


get : String -> SearchResult -> String
get key searchResult =
    Maybe.withDefault "" (Dict.get key searchResult.data)


tableConfig : Table.Config SearchResult Msg
tableConfig =
    let
        getId =
            .id >> String.fromInt
    in
    Table.config
        { toId = getId
        , toMsg = SetResultsTableState
        , columns =
            [ Table.stringColumn "Row" getId
            , Table.stringColumn "SRA Run" (get "SRR_accession")
            , Table.stringColumn "QC summary" (get "QC_summary")
            , Table.stringColumn "SRA experiment" (get "SRX_accession")
            , Table.stringColumn "SRA sample" (get "SRS_accession")
            , Table.stringColumn "SRA project" (get "SRP_accession")
            , Table.stringColumn "Sample" (get "Sample_name")
            , Table.stringColumn "Experiment" (get "GEO_series")
            ]
        }

--# This is the header of the current search result page
--# ['SRA run accession', 'QC summary alttext ', 'SRA experiment accession', 'SRA sample accession',
--# 'SRA project accession', 'Sample Name / GEO sample accession', 'GEO series accession', 'Experiment name']

viewSearchResults : Model -> Html Msg
viewSearchResults { searchResults, resultsTableState, resultsTableQuery } =
    Table.view tableConfig resultsTableState (Array.toList searchResults)
